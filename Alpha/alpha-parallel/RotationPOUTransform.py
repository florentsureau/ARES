#!/usr/bin/env python3

import numpy as np
import healpy as hp

from scipy.ndimage.filters import gaussian_filter
from math import sqrt, pi

from tqdm import tqdm

from AlphaTransformParallel import ParallelAlphaTransform as AST_Par
from AlphaTransformParallel import DTYPES, NUM_FFTW_THREADS, sharedmem_pool

# FIXME: We need to import the 'parallel' functions!
#        Later on, we can probably also remove the 'sequential' imports!
from DenoiseUtility import get_faces, set_faces, set_faces_par
from DenoiseUtility import rotate_map, rotate_map_inv


import numexpr as ne
import sharedmem as sm
import math

def cutoff_good(start, rise_length, constant_length, decrease_length):
    def actual_cutoff_good(x):
        def cutoff_help(x):
            def cutoff_help_inner(y):
                with np.errstate(divide='ignore', invalid='ignore'):
                    reciprocal = 1 / np.abs(y)
                    reciprocal[reciprocal == np.inf] = 0
                    reciprocal = np.nan_to_num(reciprocal)
                    return (y > 0) * np.exp(- reciprocal)

            temp = cutoff_help_inner(x)
            return temp / (temp + cutoff_help_inner(1 - x))

        length = rise_length + constant_length + decrease_length
        after_start = (x > start)
        after_rise = (x > start + rise_length)
        after_const = (x > start + rise_length + constant_length)
        after_decrease = (x > start + length)

        return ((after_start
                 * (1 - after_rise)
                 * cutoff_help((x - start) / rise_length))
                + (after_rise * (1 - after_const))
                + (after_const
                   * (1 - after_decrease)
                   * cutoff_help((start + length - x) / decrease_length)))

    return actual_cutoff_good

def cutoff_cos(start, rise_length, constant_length, decrease_length):
    def cutoff_func(x):
        #print('using cosine mask')
        def _cos(x):
            return -0.5*np.cos(pi*x)+0.5

        length = rise_length + constant_length + decrease_length
        after_start = (x > start)
        after_rise = (x > start + rise_length)
        after_const = (x > start + rise_length + constant_length)
        after_decrease = (x > start + length)

        return ((after_start
                 * (1 - after_rise)
                 * _cos((x - start) / rise_length))
                + (after_rise * (1 - after_const))
                + (after_const
                   * (1 - after_decrease)
                   * _cos((start + length - x) / decrease_length)))

    return cutoff_func

def cutoff_sin(start, rise_length, constant_length, decrease_length):
    def cutoff_func(x):
        #print('using sine mask')
        def _sin(x):
            return np.sin(pi*x/2)

        length = rise_length + constant_length + decrease_length
        after_start = (x > start)
        after_rise = (x > start + rise_length)
        after_const = (x > start + rise_length + constant_length)
        after_decrease = (x > start + length)

        return ((after_start
                 * (1 - after_rise)
                 * _sin((x - start) / rise_length))
                + (after_rise * (1 - after_const))
                + (after_const
                   * (1 - after_decrease)
                   * _sin((start + length - x) / decrease_length)))

    return cutoff_func

def cutoff_linear(start, rise_length, constant_length, decrease_length):
    def actual_cutoff_linear(x):
        length = rise_length + constant_length + decrease_length
        after_start = (x > start)
        after_rise = (x > start + rise_length)
        after_const = (x > start + rise_length + constant_length)
        after_decrease = (x > start + length)

        return ((after_start
                 * (1 - after_rise)
                 * ((x - start) / rise_length))
                + (after_rise * (1 - after_const))
                + (after_const
                   * (1 - after_decrease)
                   * ((start + length - x) / decrease_length)))

    return actual_cutoff_linear


class PWPOU_TransformParallel:
    def __init__(self,
                 nside,
                 alphas,
                 margin=None,
                 transition=None,
                 nested=True,
                 real=True,
                 num_cores=None,
                 num_rots=None,
                cutoff='trig'):
        self.__nside = nside
        self.__alphas = alphas
        self.__nested = nested
        self.__real = real
        self.__num_cores = num_cores
        self.__num_rots = num_rots
        self.__cutoff=cutoff
        self.__margin = margin if margin else nside // 16 - 1
        self.__transition = transition if transition else nside // 16
        self.__overlap = self.__margin + self.__transition
        self.__innermask=None
        self.__outermask=None
        assert self.__overlap <= self.__nside // 4

        # FIXME: Change this if compose_polar_faces is changed!
        self.__patch_nside = nside

        self.__trafo = AST_Par(self.__patch_nside,
                               self.__patch_nside,
                               self.__alphas,
                               real=self.__real,
                               total_cores=num_cores)
    def _create_masks(self,cutoff=None):
        vanish_length = self.__margin
        rise_length = self.__transition
        nside = self.__nside
        if cutoff is None:
            cutoff=self.__cutoff

        if cutoff is 'linear':
            inner_mask_func = cutoff_linear(-1, # no vanishing part
                                    vanish_length, # deliberate
                                    nside - 2 * (vanish_length),
                                    vanish_length) # deliberate

            outer_mask_func = cutoff_linear(vanish_length,
                                    rise_length,
                                    nside - 2 * (vanish_length + rise_length + 1),
                                    rise_length)
        elif cutoff is 'cos':
            inner_mask_func = cutoff_cos(-1, # no vanishing part
                                    vanish_length, # deliberate
                                    nside - 2 * (vanish_length),
                                    vanish_length) # deliberate

            outer_mask_func = cutoff_cos(vanish_length,
                                    rise_length,
                                    nside - 2 * (vanish_length + rise_length + 1),
                                    rise_length)
        elif cutoff is 'trig':
            inner_mask_func = cutoff_sin(-1, # no vanishing part
                                    vanish_length, # deliberate
                                    nside - 2 * (vanish_length),
                                    vanish_length) # deliberate

            outer_mask_func = cutoff_cos(vanish_length,
                                    rise_length,
                                    nside - 2 * (vanish_length + rise_length + 1),
                                    rise_length)
        else:
            inner_mask_func = cutoff_good(-1, # no vanishing part
                                    vanish_length, # deliberate
                                    nside - 2 * (vanish_length),
                                    vanish_length) # deliberate

            outer_mask_func = cutoff_good(vanish_length,
                                    rise_length,
                                    nside - 2 * (vanish_length + rise_length + 1),
                                    rise_length)

        temp_vec = inner_mask_func(np.arange(0, nside))
        inner_mask = temp_vec.reshape((nside, 1)) * temp_vec.reshape((1, nside))
        temp_vec = outer_mask_func(np.arange(0, nside))
        outer_mask = temp_vec.reshape((nside, 1)) * temp_vec.reshape((1, nside))
        return inner_mask,outer_mask

    def create_mask_old_cover(self,nside,num_rots=None):
        rotations = [[0.0, 0.0, 0.0],
                    [math.pi/4.0,0,0],
                    [0,math.pi/2.0,0],
                    [-math.pi/4.0,0,0],
                    [0,-math.pi/2.0,0],
                    [math.pi/2.,math.pi/2.0,0],
                    [math.pi/2.,-math.pi/2.0,0]]
        self.__innermask,self.__outermask=self._create_masks()
        if num_rots is None:
            NRots=1
        else:
            NRots=min((num_rots,len(rotations)))

        outer_mask_faces = np.repeat(self.__outermask.reshape((1, nside, nside)), 12, axis=0)
        assert outer_mask_faces.shape == ((12, nside, nside))
        outer_mask_map = set_faces(outer_mask_faces, nested=self.__nested)
        mask_cover = np.zeros((hp.nside2npix(nside)))
        for rot in rotations[0:NRots]:
            if [0.0] * 3 == rot:
                mask_cover += outer_mask_map
            else:
                mask_cover += rotate_map(outer_mask_map,
                                    rot[0], rot[1], rot[2],
                                    X=False, Y=True, nested=self.__nested)
        return mask_cover

    def create_mask_cover(self,nside,num_rots=None):
        if(num_rots is None):
            c = (346 * nside)//2048
            rotations = [[0.0, 0.0, 0.0],
                    [-c*pi/nside, 0.0, 0.0],
                    [c*pi/nside, 0.0, 0.0],
                    [pi/4, pi/2, 0.0],
                    [-pi/4, pi/2,0.0]]
        elif(num_rots == 1):
            rotations = [[0.0, 0.0, 0.0]]
        elif(num_rots == 4):
            c = (341 * nside)//2048
            rotations = [[0.0, 0.0, 0.0],
                [c*pi/nside, 0.0, 0.0],
                [-c*pi/nside, 0.0, 0.0],
                [0.0, pi/2, 0.0]]
        elif(num_rots== 5):
            c = (346 * nside)//2048
            rotations = [[0.0, 0.0, 0.0],
                        [c*pi/nside, 0.0, 0.0],
                        [-c*pi/nside, 0.0, 0.0],
                        [pi/4, pi/2, 0.0],
                        [-pi/4, pi/2,0.0]]
        elif(num_rots == 6):
            rotations = [[0.0, 0.0, 0.0],
                [pi/4, 0.0, 0.0],
                [pi/16, 0.0, 0.0],
                [-pi/16, 0.0, 0.0],
                [pi/4, pi/2, 0.0],
                [-pi/4, pi/2,0.0]]
        else:
            rotations = [[0.0, 0.0, 0.0],
                [pi/4, 0.0, 0.0],
                [pi/16, 0.0, 0.0],
                [-pi/16, 0.0, 0.0],
                [pi/8, 0.0, 0.0],
                [-pi/8, 0.0, 0.0],
                [pi/4, pi/2, 0.0],
                [-pi/4, pi/2,0.0]]
        self.__innermask,self.__outermask=self._create_masks()

        outer_mask_faces = np.repeat(self.__outermask.reshape((1, nside, nside)), 12, axis=0)
        assert outer_mask_faces.shape == ((12, nside, nside))
        outer_mask_map = set_faces(outer_mask_faces, nested=self.__nested)
        mask_cover = np.zeros((hp.nside2npix(nside)))
        for rot in rotations:
            if [0.0] * 3 == rot:
                mask_cover += outer_mask_map
            else:
                mask_cover += rotate_map(outer_mask_map,
                                    rot[0], rot[1], rot[2],
                                    X=False, Y=True, nested=self.__nested)
        return mask_cover


    def _denoisePOU(self, hp_im,noise_sigma=None,multipliers=None):
        nside=self.__nside
        if(self.__num_rots is None):
            c = (346 * nside)//2048
            rotations = [[0.0, 0.0, 0.0],
                    [-c*pi/nside, 0.0, 0.0],
                    [c*pi/nside, 0.0, 0.0],
                    [pi/4, pi/2, 0.0],
                    [-pi/4, pi/2,0.0]]
        elif(self.__num_rots == 1):
            rotations = [[0.0, 0.0, 0.0]]
        elif(self.__num_rots == 3): #ORIGINAL CYCLE SPAN
            rotations = [[0.0, 0.0, 0.0],
                [math.pi/4.0,0,0],
                [0,math.pi/2.0,0],
                [-math.pi/4.0,0,0],
                [0,-math.pi/2.0,0],
                [math.pi/2.,math.pi/2.0,0],
                [math.pi/2.,-math.pi/2.0,0]]
        elif(self.__num_rots == 4):
            c = (341 * nside)//2048
            rotations = [[0.0, 0.0, 0.0],
                [c*pi/nside, 0.0, 0.0],
                [-c*pi/nside, 0.0, 0.0],
                [0.0, pi/2, 0.0]]
        elif(self.__num_rots == 5):
            c = (346 * nside)//2048
            rotations = [[0.0, 0.0, 0.0],
                        [c*pi/nside, 0.0, 0.0],
                        [-c*pi/nside, 0.0, 0.0],
                        [pi/4, pi/2, 0.0],
                        [-pi/4, pi/2,0.0]]
        elif(self.__num_rots == 6):
            rotations = [[0.0, 0.0, 0.0],
                [pi/4, 0.0, 0.0],
                [pi/16, 0.0, 0.0],
                [-pi/16, 0.0, 0.0],
                [pi/4, pi/2, 0.0],
                [-pi/4, pi/2,0.0]]
        else:
            rotations = [[0.0, 0.0, 0.0],
                [pi/4, 0.0, 0.0],
                [pi/16, 0.0, 0.0],
                [-pi/16, 0.0, 0.0],
                [pi/8, 0.0, 0.0],
                [-pi/8, 0.0, 0.0],
                [pi/4, pi/2, 0.0],
                [-pi/4, pi/2,0.0]]
        self.__innermask,self.__outermask=self._create_masks()

        result = np.zeros_like(hp_im)
        with tqdm(total=len(rotations), desc='Denoising faces', ncols=80, unit='rot') as pbar:
            for rot in rotations:
                if [0.0] * 3 == rot:
                    print('----> rotation is the identity')
                    coeffs = self.transform_generator(hp_im)
                    faces=self.denoise_transform(coeffs,multipliers=multipliers,noise_sigma=noise_sigma)
                    result = set_faces(faces, nested=self.__nested)
                else:
                    print('----> rotation is ',rot)
                    rotmap = rotate_map_inv(hp_im,rot[0], rot[1], rot[2],
                                    X=False, Y=True, nested=self.__nested)
                    coeffs = self.transform_generator(rotmap)
                    faces=self.denoise_transform(coeffs,multipliers=multipliers,noise_sigma=noise_sigma)
                    result += rotate_map(set_faces(faces, nested=self.__nested),
                                         rot[0], rot[1], rot[2],
                                         X=False, Y=True, nested=self.__nested)
                pbar.update()

        self.__trafo = None # clear the memory to make space for performing the rotations
        outer_mask_faces = np.repeat(self.__outermask.reshape((1, nside, nside)), 12, axis=0)
        assert outer_mask_faces.shape == ((12, nside, nside))
        outer_mask_map = set_faces(outer_mask_faces, nested=self.__nested)
        mask_sum = np.zeros_like(hp_im)
        for rot in rotations:
            if [0.0] * 3 == rot:
                mask_sum += outer_mask_map
            else:
                mask_sum += rotate_map(outer_mask_map,
                                    rot[0], rot[1], rot[2],
                                    X=False, Y=True, nested=self.__nested)

        return result, mask_sum

    def transform_generator(self, hp_im, do_norm=True):
        #ITERATOR ON THE TRANSFORM
        print("TRANS GEN")
        assert hp.get_nside(hp_im) == self.__nside
        faces = get_faces(hp_im, nested=self.__nested)
        for kf in range(12):
            yield self.__trafo.transform(faces[kf], do_norm=do_norm)

    def denoise_transform(self,coeffs_gen,noise_sigma=None,multipliers=None):
        #TAKE SEQUENTIALLY ALL COEFFS BELONGING TO A FACE, THR THEM, RECONSTRUCT
        nside=self.__nside
        faces = np.zeros((12, nside, nside))
        for i, coeffs in zip(range(12), coeffs_gen):
            print("PROCESS eq face",i)
            faces[i] = self.__outermask*self.inv_thresh_scales(coeffs,noise_sigma=noise_sigma,
                                                     multipliers=multipliers)
        return faces

    def inv_thresh_scales(self,coeffs, noise_sigma=None,multipliers=None):
        #INV THRESHOLD THE SCALES FOR A FACE
        if noise_sigma is None:
            print('ABORT: a sigma value must be specified')
            return hp_im
        assert noise_sigma >= 0.0

        if multipliers is None:
            multipliers = [3] * trafo.num_scales + [4]
        width = self.__trafo.width
        height = self.__trafo.height
        indices=self.__trafo.indices
        for kl,ksc in enumerate(indices):
            if(kl==0):
                scales=[0]
            else:
                scales.append(ksc[0]+1)
        thresh_lvls = [noise_sigma*multipliers[sc] for sc in scales]
        thresh_coeff=sm.empty_like(coeffs,dtype=DTYPES[1])
        with sharedmem_pool(self.__num_cores) as pool:
            def work(i):
                coeff=coeffs[i]
                thresh=thresh_lvls[i]
                ne.evaluate('coeff*(real(abs(coeff))>=thresh)',
                                                out=thresh_coeff[i])
            pool.map(work,range(len(indices)))

        thresh_coeff = np.array(thresh_coeff)
        recon =self.__trafo.inverse_transform(thresh_coeff, real=True, do_norm=True)
        return recon

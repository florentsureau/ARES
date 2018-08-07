#!/usr/bin/env python3

import numpy as np
import healpy as hp

from tqdm import tqdm

from AlphaTransformParallel import ParallelAlphaTransform as AST_Par
from AlphaTransformParallel import DTYPES, NUM_FFTW_THREADS, sharedmem_pool

from DenoiseUtility import get_faces, set_faces, set_faces_par
from DenoiseUtility import compose_polar_faces, compose_equatorial_faces
from DenoiseUtility import reassemble_faces, reassemble_faces_par

import numexpr as ne
import sharedmem as sm
import math

class PW_TransformParallel:
    def __init__(self,
                 nside,
                 alphas,
                 margin=None,
                 transition=None,
                 nested=True,
                 real=True,
                 num_cores=None):
        self.__nside = nside
        self.__alphas = alphas
        self.__nested = nested
        self.__real = real
        self.__num_cores = num_cores

        self.__margin = margin if margin else nside // 16
        self.__transition = transition if transition else nside // 16
        self.__overlap = self.__margin + self.__transition
        assert self.__overlap <= self.__nside // 4

        # FIXME: Change this if compose_polar_faces is changed!
        self.__patch_nside = 2 * (3*nside//4 + self.__overlap)

        self.__trafo = AST_Par(self.__patch_nside,
                               self.__patch_nside,
                               self.__alphas,
                               real=self.__real,
                               total_cores=num_cores)

    def _face_generator(self, hp_im):
        assert hp.get_nside(hp_im) == self.__nside
        faces = get_faces(hp_im, nested=self.__nested)

        equat_composite_faces = compose_equatorial_faces(faces,
                                                         self.__overlap,
                                                         self.__num_cores)
        assert equat_composite_faces[0].shape[0] == self.__patch_nside
        assert equat_composite_faces[0].shape[1] == self.__patch_nside

        yield from equat_composite_faces
        del equat_composite_faces

        polar_composite_faces = compose_polar_faces(faces,
                                                    self.__overlap,
                                                    self.__num_cores)

        assert polar_composite_faces[0].shape[0] == self.__patch_nside
        assert polar_composite_faces[0].shape[1] == self.__patch_nside

        yield from polar_composite_faces


    def transform_generator(self, hp_im, do_norm=True):
        for face in self._face_generator(hp_im):
            yield self.__trafo.transform(face, do_norm=do_norm)

    def max_coeff(self, hp_im, do_norm=True):
        maxima = [0 for j in range(6)]

        for j, face in enumerate(self._face_generator(hp_im)):
            face_fourier = self.__trafo.do_fourier(face)

            numexpr_flag = NUM_FFTW_THREADS != 1
            with sharedmem_pool(self.__num_cores, numexpr=numexpr_flag) as pool:
                def work(i):
                    coeff = self.__trafo.transform_fourier(face_fourier,
                                                           i,
                                                           do_norm=do_norm)
                    return np.max(np.abs(coeff))

                maxima[j] = np.max(pool.map(work,
                                            range(self.__trafo.redundancy)))

        return np.max(maxima)

    def inverse_transform(self, coeffs_gen, real=True, do_norm=True):
        eq_faces = [None, None, None, None]

        for i, coeffs in zip(range(4), coeffs_gen):
            eq_faces[i] = self.__trafo.inverse_transform(coeffs,
                                                         real=real,
                                                         do_norm=do_norm)

        polar_faces = [None, None]

        for i, coeffs in zip(range(2), coeffs_gen):
            polar_faces[i] = self.__trafo.inverse_transform(coeffs,
                                                            real=real,
                                                            do_norm=do_norm)

        faces = reassemble_faces_par(self.__nside,
                                     polar_faces,
                                     eq_faces,
                                     self.__margin,
                                     self.__transition,
                                     self.__num_cores)

        return set_faces(faces, nested=self.__nested)

    def inverse_thresholded_transform(self,
                                      coeffs_gen,
                                      threshold,
                                      real=True,
                                      do_norm=True):
        eq_faces = [None, None, None, None]

        for i, coeffs in zip(range(4), coeffs_gen):
            dn = do_norm
            eq_faces[i] = self.__trafo.inverse_thresholded_transform(coeffs,
                                                                     threshold,
                                                                     real=real,
                                                                     do_norm=dn)

        pol_fcs = [None, None]

        for i, coeffs in zip(range(2), coeffs_gen):
            dn = do_norm
            pol_fcs[i] = self.__trafo.inverse_thresholded_transform(coeffs,
                                                                    threshold,
                                                                    real=real,
                                                                    do_norm=dn)

        faces = reassemble_faces_par(self.__nside,
                                     pol_fcs,
                                     eq_faces,
                                     self.__margin,
                                     self.__transition,
                                     self.__num_cores)

        return set_faces(faces, nested=self.__nested)

    def denoise_transform(self,coeffs_gen,noise_sigma=None,multipliers=None):
        eq_faces = [None, None, None, None]

        for i, coeffs in zip(range(4), coeffs_gen):
            print("PROCESS eq face",i)
            eq_faces[i] = self.inv_thresh_scales(coeffs,noise_sigma=noise_sigma,
                                                     multipliers=multipliers)

        pol_fcs = [None, None]

        for i, coeffs in zip(range(2), coeffs_gen):
            print("PROCESS polar face",i)
            pol_fcs[i] = self.inv_thresh_scales(coeffs,noise_sigma=noise_sigma,
                                                     multipliers=multipliers)

        faces = reassemble_faces_par(self.__nside,
                                     pol_fcs,
                                     eq_faces,
                                     self.__margin,
                                     self.__transition,
                                     self.__num_cores)

        return set_faces(faces, nested=self.__nested)

    def inv_thresh_scales(self,coeffs, noise_sigma=None,multipliers=None):
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


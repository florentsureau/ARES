#!/usr/bin/env python3
# pylint: disable=fixme

import math
import os
import numpy as np
import numexpr as ne
import pyfftw

from time import perf_counter
from tqdm import tqdm

import sharedmem as sm
# import multiprocessing

import MotherShearlets as MS
from fourier_util import fft2, ifft2, my_fft_shift, my_ifft_shift

# FIXME: NOTE: I added 'numexpr=False' to some calls of 'sharedmem_pool'.
#              In these calls, numexpr is not used, but usually pyfftw is.
#              Thus, if one really uses threaded pyfftw,
#              then this might need to be changed!

# DTYPES = ('f4', 'complex64') # for single precision
DTYPES = ('f8', 'complex128') # for double precision

def div0(a, b):
    r"""
    This (numpy-vectorized) function computes a/b, ignoring divisions by 0.

    Example:
        div0( [-1, 0, 1], 0 ) -> [0, 0, 0]
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[~ np.isfinite(c)] = 0  # -inf inf NaN
    return c

def good_process_number(total_cores, numexpr=True):
    if numexpr:
        thread_num = ne.detect_number_of_threads()
    else:
        thread_num = 1
    return max(1, total_cores // thread_num)


def calculate_bounds(width, height):
    epsilon_width = 1 - (width % 2)
    # Note: the result is always an integer,
    # but python does not know this, so we use // division
    n = (width - 1 - epsilon_width) // 2
    assert width == 2 * n + 1 + epsilon_width

    epsilon_height = 1 - (height % 2)
    m = (height - 1 - epsilon_height) // 2
    assert height == 2 * m + 1 + epsilon_height

    return (-n - epsilon_width, n, -m - epsilon_height, m)


def build_rescaled_functions(x_min,
                             x_max,
                             y_min,
                             y_max,
                             mother_shearlet,
                             alphas):
    num_scales = len(alphas)
    # calculate (horizontal and vertical) scaling factors
    n_max = - x_min
    m_max = - y_min

    scale_function = mother_shearlet.scale_function
    scale_fun_lower_bound = scale_function.large_support[0]
    scale_fun_upper_bound = scale_function.large_support[1]

    max_scale = 2 ** (-(num_scales - 1))

    R = min(n_max, m_max)
    # this is the size of the (quadratic(!)) low pass region
    a = max_scale * scale_fun_lower_bound / scale_fun_upper_bound * R
    assert a >= 4, ("The given number of scales is too large "
                    "for the given dimensions!")

    scale_nominator = ((scale_fun_upper_bound - scale_fun_lower_bound) /
                       max_scale)
    hor_scale = scale_nominator / (n_max - a / max_scale)
    vert_scale = scale_nominator / (m_max - a / max_scale)
    hor_shift = hor_scale * a - scale_fun_lower_bound
    vert_shift = vert_scale * a - scale_fun_lower_bound

    horizontal_scale_fct = MS.scale(MS.translate(scale_function, hor_shift),
                                    1 / hor_scale)
    vertical_scale_fct = MS.scale(MS.translate(scale_function, vert_shift),
                                  1 / vert_scale)

    # rescale the low-pass functions
    orig_low_pass_fct = mother_shearlet.low_pass_function

    hori_scale_low_pass = (a / orig_low_pass_fct.large_support[1])
    horizontal_low_pass_fct = MS.scale(orig_low_pass_fct, hori_scale_low_pass)

    vert_scale_low_pass = (a / orig_low_pass_fct.large_support[1])
    vert_low_pass_fct = MS.scale(orig_low_pass_fct, vert_scale_low_pass)

    return (horizontal_scale_fct,
            vertical_scale_fct,
            horizontal_low_pass_fct,
            vert_low_pass_fct)


AFFINITY_FLAG = False

def sharedmem_pool(total_cores, numexpr=True):
    # see https://stackoverflow.com/questions/15639779
    global AFFINITY_FLAG
    if not AFFINITY_FLAG:
        AFFINITY_FLAG = True
        os.system("taskset -p 0xfff %d" % os.getpid())
    if total_cores is None:
        total_cores = sm.cpu_count()
    return sm.MapReduce(np=good_process_number(total_cores, numexpr))
    # return sm.MapReduceByThread(np=good_process_number())


def xy_grids(x_min, x_max, y_min, y_max):
    x_values = np.mgrid[x_min:x_max + 1]
    y_values = np.mgrid[y_max:y_min - 1:-1].reshape((y_max - y_min + 1, 1))

    return (x_values, y_values)


def calc_grids(width, height, periodization):
    (x_min, x_max, y_min, y_max) = calculate_bounds(width, height)
    (values_x, values_y) = xy_grids(x_min, x_max, y_min, y_max)

    if not periodization:
        x_values = [values_x]
        y_values = [values_y]

        hor_quotient_grids = [div0(values_y, values_x)]
        vert_quotient_grids = [div0(values_x, values_y)]
    else:
        xy_values = [xy_grids(x_min + i * width, x_max + i * width,
                              y_min + j * height, y_max + j * height)
                     for i in range(-1, 2)
                     for j in range(-1, 2)]

        x_values = [value[0] for value in xy_values]
        y_values = [value[1] for value in xy_values]
        hor_quotient_grids = [div0(y_value, x_value)
                              for x_value, y_value in zip(x_values, y_values)]
        vert_quotient_grids = [div0(x_value, y_value)
                               for x_value, y_value in zip(x_values, y_values)]

    return (values_x, values_y,
            x_values, y_values,
            hor_quotient_grids, vert_quotient_grids)


def calculate_single_spect(width,
                           height,
                           x_values,
                           y_values,
                           hor_quo_grids,
                           ver_quo_grids,
                           index,
                           alphas,
                           direction_function,
                           hor_scale_fct,
                           vert_scale_fct,
                           real):
    # 'unnormalized' can be used to overwrite the "parseval" normalization.
    # This is important, since the parseval normalization requires the
    # dual_frame_weight, which can only be computed once all "unnormalized"
    # spectrograms are known!
    (j, k, cone) = index

    cur_hor_scale_fct = MS.scale(hor_scale_fct, 2**j)
    cur_vertical_scale_fct = MS.scale(vert_scale_fct, 2**j)
    dir_fct = MS.scale(MS.translate(direction_function, k),
                       2**((alphas[j] - 1) * j))

    # cur_spect = np.zeros((height, width))
    cur_spect = np.zeros((height, width), dtype=DTYPES[0])

    for (x_value, y_value,
         hor_quo_grid, vert_quo_grid) in zip(x_values, y_values,
                                             hor_quo_grids, ver_quo_grids):
        # fac2 = np.zeros(hor_quo_grid.shape)
        fac2 = np.zeros(hor_quo_grid.shape, dtype=DTYPES[0])

        if real:
            if cone == 'h':
                first_fact1 = cur_hor_scale_fct.call(x_value)
                first_fact2 = cur_hor_scale_fct.call(-x_value)
                fac1 = first_fact1 + first_fact2
                fac2[:, fac1 > 0] = dir_fct.call(hor_quo_grid[:, fac1 > 0])
            elif cone == 'v':
                first_fact1 = cur_vertical_scale_fct.call(y_value)
                first_fact2 = cur_vertical_scale_fct.call(-y_value)
                fac1 = first_fact1 + first_fact2
                indices = (fac1 > 0).ravel()
                fac2[indices, :] = dir_fct.call(vert_quo_grid[indices, :])

            ne.evaluate('cur_spect + fac1 * fac2', out=cur_spect)

        else:
            if cone == 'r':
                fac1 = cur_hor_scale_fct.call(x_value)
                fac2[:, fac1 > 0] = dir_fct.call(hor_quo_grid[:, fac1 > 0])
            elif cone == 't':
                fac1 = cur_vertical_scale_fct.call(y_value)
                indices = (fac1 > 0).ravel()
                fac2[indices, :] = dir_fct.call(vert_quo_grid[indices, :])
            elif cone == 'l':
                fac1 = cur_hor_scale_fct.call(-x_value)
                fac2[:, fac1 > 0] = dir_fct.call(hor_quo_grid[:, fac1 > 0])
            elif cone == 'b':
                fac1 = cur_vertical_scale_fct.call(-y_value)
                assert fac1 is not None  # silence pyflakes
                indices = (fac1 > 0).ravel()
                fac2[indices, :] = dir_fct.call(vert_quo_grid[indices, :])
                assert fac2 is not None  # silence pyflakes

            ne.evaluate('cur_spect + fac1 * fac2', out=cur_spect)

    return cur_spect


def compute_dual_weight_and_norms(width, height, spectrograms, total_cores):
    # currently, this is not parallelized!
    # FIXME: IT COULD BE A GOOD IDEA TO CHANGE THIS!
    # dual_frame_weight = np.zeros((height, width))
    dual_frame_weight = np.zeros((height, width), dtype=DTYPES[0])
    for spect in spectrograms:
        dual_frame_weight += np.square(spect)

    # fourier_norms = sm.empty(len(spectrograms))
    # space_norms = sm.empty(len(spectrograms))
    fourier_norms = sm.empty(len(spectrograms), dtype=DTYPES[0])
    space_norms = sm.empty(len(spectrograms), dtype=DTYPES[0])

    with sharedmem_pool(total_cores, numexpr=False) as pool:
        def work(i):
            spect_norm = np.linalg.norm(spectrograms[i])
            fourier_norms[i] = spect_norm
            space_norms[i] = spect_norm / math.sqrt(spectrograms[i].size)

        pool.map(work, list(range(len(spectrograms))))

    return (dual_frame_weight, list(fourier_norms), list(space_norms))


# This function works in parallel, and either creates a new shared memory
# array into which the spectrograms are saved.
#
# Return value: A tuple (spects, dual_frame_weights), where the shared memory
# array 'spects' contains the computed spectrograms, and 'dual_frame_weight'
# contains the computed dual frame weights.
#
# This method assumes that:
# * we want fully sampled spectrograms,
# * we use the standard mother shearlet
def calculate_full_spects(width,
                          height,
                          alphas,
                          real,
                          parseval,
                          periodization,
                          total_cores):
    indices = ParallelAlphaTransform.calculate_indices(real, alphas)
    # Note: the spects are real (at least with the Haeuser mother shearlet)
    # Note: the original implementation also uses double precision
    #       (i.e. float64 alias 'f8')
    # Note: Idea for saving memory in the future:
    #       USE SINGLE PRECISION INSTEAD OF DOUBLE!
    # spects = sm.empty((len(indices), height, width))
    spects = sm.empty((len(indices), height, width), dtype=DTYPES[0])

    (x_min, x_max, y_min, y_max) = calculate_bounds(width, height)
    num_scales = len(alphas)

    mother_shearlet = MS.HaeuserMotherShearlet
    (hor_scale_fct,
     vert_scale_fct,
     hor_low_pass,
     vert_low_pass) = build_rescaled_functions(x_min, x_max,
                                               y_min, y_max,
                                               mother_shearlet,
                                               alphas)

    direction_function = mother_shearlet.direction_function

    (values_x, values_y,
     x_values, y_values,
     hor_quo_grids, vert_quo_grids) = calc_grids(width, height, periodization)

    spects[0] = hor_low_pass.call(values_x) * vert_low_pass.call(values_y)
    del values_x
    del values_y

    # spect_start = perf_counter()
    enumerated_indices = list(enumerate(indices))[1:]
    with sharedmem_pool(total_cores) as pool:
        def work(par):
            (i, index) = par
            spects[i] = calculate_single_spect(width, height,
                                               x_values, y_values,
                                               hor_quo_grids, vert_quo_grids,
                                               index,
                                               alphas,
                                               direction_function,
                                               hor_scale_fct,
                                               vert_scale_fct,
                                               real)

        pool.map(work, enumerated_indices)

    # spect_end = perf_counter()
    # print(spect_end - spect_start)

    (dual_frame_weight,
     fourier_norms,
     space_norms) = compute_dual_weight_and_norms(width, height,
                                                  spects, total_cores)

    if parseval:
        with sharedmem_pool(total_cores) as pool:
            def work(i):
                spec = spects[i]
                dual_weight = dual_frame_weight
                ne.evaluate('spec / sqrt(dual_weight)', out=spects[i])
                # spects[i] /= np.sqrt(dual_frame_weight)

            pool.map(work, list(range(len(indices))))

    return (np.array(spects), dual_frame_weight, fourier_norms, space_norms)

# FIXME: In general, I recommend setting this to ne.detect_number_of_threads()
#        But on my computer, fftw() gets stuck if I use more than one thread...
NUM_FFTW_THREADS = 1


def make_fft(width, height):
    # threads = ne.detect_number_of_cores()
    # threads = ne.detect_number_of_threads()
    threads = NUM_FFTW_THREADS
    # print("Building pyfftw object with threads = {0}".format(threads))
    # inp = pyfftw.empty_aligned((height, width), dtype='complex128')
    # out = pyfftw.empty_aligned((height, width), dtype='complex128')
    inp = pyfftw.empty_aligned((height, width), dtype=DTYPES[1])
    out = pyfftw.empty_aligned((height, width), dtype=DTYPES[1])
    return pyfftw.FFTW(inp, out, axes=(0, 1), threads=threads)

def make_ifft(width, height):
    # threads = ne.detect_number_of_cores()
    threads = NUM_FFTW_THREADS
    # print("Building pyfftw object with threads = {0}".format(threads))
    # inp = pyfftw.empty_aligned((height, width), dtype='complex128')
    # out = pyfftw.empty_aligned((height, width), dtype='complex128')
    inp = pyfftw.empty_aligned((height, width), dtype=DTYPES[1])
    out = pyfftw.empty_aligned((height, width), dtype=DTYPES[1])
    return pyfftw.FFTW(inp,
                       out,
                       axes=(0, 1),
                       direction='FFTW_BACKWARD',
                       threads=threads)

class ParallelAlphaTransform:
    def __init__(self,
                 width,
                 height,
                 alphas,
                 *,
                 # if 'None' is used, then actually sm.cpu_count() is used!
                 total_cores=None,
                 real=True,
                 parseval=False,
                 periodization=True):
        assert all(0 <= alpha <= 1 for alpha in alphas), ("alphas must be "
                                                          "between 0 and 1.")

        # Set the relevant parameters
        self._width = width
        self._height = height
        self._alphas = alphas
        self._total_cores = total_cores
        self._real = real
        self._parseval = parseval
        self._periodization = periodization
        self._indices = ParallelAlphaTransform.calculate_indices(real, alphas)

        (spectrograms,
         dual_frame_weight,
         fourier_norms,
         space_norms) = calculate_full_spects(width,
                                              height,
                                              alphas,
                                              real,
                                              parseval,
                                              periodization,
                                              total_cores)

        self._spectrograms = spectrograms
        self._dual_frame_weight = dual_frame_weight
        self._fourier_norms = fourier_norms
        self._space_norms = space_norms
        self._frame_bounds = (np.min(dual_frame_weight),
                              np.max(dual_frame_weight))

        self._setup_fft()

    def _setup_fft(self):
        # threads = ne.detect_number_of_cores()
        # height = self._height
        # width = self._width

        # inp = pyfftw.empty_aligned((height, width), dtype=DTYPES[1])
        # out = pyfftw.empty_aligned((height, width), dtype=DTYPES[1])
        # fft_setup_start = perf_counter()
        # self.__fftw = pyfftw.FFTW(inp, out, axes=(0, 1), threads=threads)
        # constructions of the FFTW don't need to compute a new 'wisdom'
        # temp = pyfftw.FFTW(inp, out, axes=(0, 1), threads=threads)
        # fft_setup_end = perf_counter()
        # print("Time for setup of FFT: {0} seconds.".format(fft_setup_end - fft_setup_start))
        # del temp

        # with sharedmem_pool() as pool:
        #     def work(i):
        #         fft_setup_start = perf_counter()
        #         process_local_inp = pyfftw.empty_aligned((height, width),
        #                                                  dtype=DTYPES[1])
        #         process_local_out = pyfftw.empty_aligned((height, width),
        #                                                  dtype=DTYPES[1])
        #         process_local_fftw = pyfftw.FFTW(process_local_inp,
        #                                          process_local_out,
        #                                          axes=(0, 1),
        #                                          threads=threads)
        #         fft_setup_end = perf_counter()
        #         print_str = "Time for FFT-setup on process {0}: {1} seconds."

        #         with pool.critical:
        #             print(print_str.format(i, fft_setup_end - fft_setup_start))

        #     pool.map(work, [0,1])

        # inp = pyfftw.empty_aligned((height, width), dtype=DTYPES[1])
        # out = pyfftw.empty_aligned((height, width), dtype='complex128')
        # fft_setup_start = perf_counter()
        # temp = pyfftw.FFTW(inp,
        # self.__ifftw = pyfftw.FFTW(inp,
        #                            out,
        #                            axes=(0, 1),
        #                            direction='FFTW_BACKWARD',
        #                            threads=threads)
        # del temp
        # fft_setup_end = perf_counter()
        # print("Time for setup of iFFT: {0} seconds.".format(fft_setup_end - fft_setup_start))
        self.__fftw = make_fft(self.width, self.height)
        self.__ifftw = make_ifft(self.width, self.height)
        self.__normalization = math.sqrt(self.width * self.height)

    @classmethod
    def calculate_indices(cls, real, alphas):
        if real:
            cones = ['h', 'v']
        else:
            cones = ['r', 't', 'l', 'b']

        indices = [-1]  # -1 is the low-pass filter

        for j, alpha in enumerate(alphas):
            for cone in cones:
                # k = int(round(2**((1 - alpha) * j)))
                k = math.ceil(2**((1 - alpha) * j))

                # to get a nice visual appearance, we change
                # the ordering of the shears on certain cones.
                shear_range = range(-k, k + 1)
                if cone in ['t', 'b', 'v']:
                    shear_range = range(k, -k - 1, -1)

                for shear in shear_range:
                    indices.append((j, shear, cone))
        return indices

    # def _fft(self, inp):
    def _fft(self, inp, fftw_obj=None):
        if fftw_obj is None:
            fftw_obj = self.__fftw
        # return self.__fftw(inp) / self.__normalization
        return fftw_obj(inp) / self.__normalization

    # def _ifft(self, inp):
    def _ifft(self, inp, ifftw_obj=None):
        if ifftw_obj is None:
            ifftw_obj = self.__ifftw
        # print("Starting '_ifft'...")
        # return self.__ifftw(inp, normalise_idft=False) / self.__normalization
        return ifftw_obj(inp, normalise_idft=False) / self.__normalization
        # print("'_ifft' finished!")

    @property
    def spectrograms(self):
        return self._spectrograms

    @property
    def shearlets(self):
        # im = np.zeros((self.height, self.width))
        im = np.zeros((self.height, self.width), dtype=DTYPES[0])
        im[self.height // 2, self.width // 2] = 1
        if self.is_real:
            return np.real(self.transform(im, do_norm=False))
        else:
            return self.transform(im, do_norm=False)

    @classmethod
    def num_shears(cls, alphas, j, real=False):
        if j == -1:
            return 1
        # shears_per_cone = 2 * int(round(2**((1 - alphas[j]) * j))) + 1
        shears_per_cone = 2 * math.ceil(2**((1 - alphas[j]) * j)) + 1
        if not real:
            return 4 * shears_per_cone
        else:
            return 2 * shears_per_cone

    def get_frame_bounds(self, do_norm=False):
        if not do_norm:
            if self.is_parseval:
                return (1, 1)
            else:
                return self._frame_bounds
        else:
            # fourier_sum = np.zeros((self.height, self.width))
            fourier_sum = np.zeros((self.height, self.width), dtype=DTYPES[0])
            for spect, norm in zip(self.spectrograms, self.space_norms):
                fourier_sum += np.square(spect / norm)
            return (np.min(fourier_sum), np.max(fourier_sum))

    @property
    def indices(self):
        return self._indices

    def scale_slice(self, scale):
        assert -1 <= scale < self.num_scales, ("The given scale is outside "
                                               "the valid range")
        if scale == -1:
            return slice(0, 1)
        lower_bound = next(i
                           for i, index in enumerate(self.indices)
                           if index != -1 and index[0] == scale)
        upper_bound = (len(self.indices) if scale + 1 == self.num_scales
                       else next(i
                                 for i, index in enumerate(self.indices)
                                 if index != -1 and index[0] == scale + 1))
        return slice(lower_bound, upper_bound)

    @property
    def periodization(self):
        return self._periodization

    @property
    def num_scales(self):
        return len(self._alphas)

    @property
    def width(self):
        return self._width

    @property
    def height(self):
        return self._height

    @property
    def is_parseval(self):
        return self._parseval

    @property
    def is_real(self):
        return self._real

    @property
    def total_cores(self):
        return self._total_cores

    @total_cores.setter
    def total_cores(self, cores):
        self._total_cores = cores

    @property
    def fourier_norms(self):
        return self._fourier_norms

    @property
    def space_norms(self):
        return self._space_norms

    @property
    def redundancy(self):
        return len(self.indices)

    def alpha(self, scale):
        return self._alphas[scale]

    def do_fourier(self, im):
        return my_fft_shift(self._fft(im))

    def _check_dimensions(self, im):
        assert (self.width == im.shape[1] and
                self.height == im.shape[0]), ("Dimensions of the "
                                              "image ({0}) do not match those "
                                              "of the transform ({1}).").format(
                                                  image.shape,
                                                  (self.height, self.width))

    # Idea of this method: Compute only the part of the transform corresponding
    #                      to self.transform(...)[i].
    #                      To avoid computing FFT of the input several times,
    #                      this method uses the Fourier-transformed version of
    #                      the input image as input.
    #                      This FFT can be computed using self.do_fourier(...).
    def transform_fourier(self, im_fourier, i, do_norm=True):
        self._check_dimensions(im_fourier)

        if do_norm:
            return (self._ifft(my_ifft_shift(im_fourier *
                                             self._spectrograms[i]))
                    / self._space_norms[i])
        else:
            return self._ifft(my_ifft_shift(im_fourier * self._spectrograms[i]))

    def transform(self, image, do_norm=True, total_cores=None):
        self._check_dimensions(image)

        if total_cores is None:
            total_cores = self.total_cores

        # im_fourier = my_fft_shift(self._fft(image))
        im_fourier = self.do_fourier(image)

        # trafo = sm.empty((len(self._spectrograms), self.height, self.width),
        #                  dtype='complex128')
        trafo = sm.empty((len(self._spectrograms), self.height, self.width),
                         dtype=DTYPES[1])
        # sm.set_debug(True)
        numexpr_flag = NUM_FFTW_THREADS != 1
        with sharedmem_pool(total_cores, numexpr=numexpr_flag) as pool:
            def work(i):
                trafo[i] = self.transform_fourier(im_fourier, i, do_norm)
                # if do_norm:
                #     # trafo[i] = (self._ifft(my_ifft_shift(im_fourier *
                #     #                                      self._spectrograms[i]),
                #     #                        local_ifft)
                #     #             / self._space_norms[i])
                #     trafo[i] = (self._ifft(my_ifft_shift(im_fourier *
                #                                          self._spectrograms[i]))
                #                 / self._space_norms[i])
                # else:
                #     # trafo[i] = self._ifft(my_ifft_shift(im_fourier *
                #     #                                     self._spectrograms[i]),
                #     #                       local_ifft)
                #     trafo[i] = self._ifft(my_ifft_shift(im_fourier *
                #                                         self._spectrograms[i]))
            pool.map(work, list(range(len(self._spectrograms))))

        return np.array(trafo)

    def add_adjoint_part(self, coeff, i, result, do_norm=True, pool=None):
        if do_norm:
            temp = ((self._spectrograms[i] / self._space_norms[i]) *
                    my_fft_shift(self._fft(coeff)))
        else:
            temp = (self._spectrograms[i] * my_fft_shift(self._fft(coeff)))

        if pool is None:
            result += temp
        else:
            with pool.critical:
                result[:] += temp

    def adjoint_transform(self, coeffs, do_norm=True, total_cores=None):
        if total_cores is None:
            total_cores = self.total_cores

        # result = sm.full((self.height, self.width), 0, dtype='complex128')
        result = sm.full((self.height, self.width), 0, dtype=DTYPES[1])

        numexpr_flag = NUM_FFTW_THREADS != 1
        with sharedmem_pool(total_cores, numexpr=numexpr_flag) as pool:
            def work(i):
                # local_fft = make_fft(self.width, self.height)
                self.add_adjoint_part(coeffs[i], i, result, do_norm, pool)
                # if do_norm:
                #     # temp = ((spects[i] / self._space_norms[i]) *
                #     #         my_fft_shift(self._fft(coeffs[i], local_fft)))
                #     temp = ((self._spectrograms[i] / self._space_norms[i]) *
                #             my_fft_shift(self._fft(coeffs[i])))
                # else:
                #     # temp = (self._spectrograms[i] *
                #     #         my_fft_shift(self._fft(coeffs[i], local_fft)))
                #     temp = (self._spectrograms[i] *
                #             my_fft_shift(self._fft(coeffs[i])))

                # with pool.critical:
                #     result[:] += temp

            pool.map(work, list(range(len(self._spectrograms))))

        # return ifft2(my_ifft_shift(result))
        return self._ifft(my_ifft_shift(result))

    @property
    def dual_frame_weight(self):
        return self._dual_frame_weight

    # This method adds the partial contribution of the i-th part of the
    # coefficients to 'result'.
    # If 'pool' is not None, but a sharedmem.MapReduce object, then the final
    # addition is guarded using a critical section.
    # If 'pool' is None, then 'result' should be an ordinary numpy array,
    # not a sharedmem.*** object!
    def add_inverse_part(self, coeff, i, result, do_norm=True, pool=None):
        coeff_fourier = my_fft_shift(self._fft(coeff))
        assert coeff_fourier is not None  # silence pyflakes
        spec = self._spectrograms[i]
        if self.is_parseval:
            dual_w = np.ones_like(self.dual_frame_weight)
        else:
            dual_w = self.dual_frame_weight

        # temp = np.zeros((self.height, self.width), dtype='complex128')
        temp = np.zeros((self.height, self.width), dtype=DTYPES[1])

        if do_norm:
            norm = self._space_norms[i]
            ne.evaluate('spec * norm * coeff_fourier / dual_w', out=temp)
        else:
            ne.evaluate('spec * coeff_fourier / dual_w', out=temp)

        if pool is None:
            result += temp
        else:
            with pool.critical:
                result[:] += temp

    def inverse_thresholded_transform(self,
                                      coeffs,
                                      threshold,
                                      real=False,
                                      do_norm=True,
                                      total_cores=None):
        thresh_coeffs = sm.empty_like(coeffs, dtype=DTYPES[1])

        if total_cores is None:
            total_cores = self.total_cores

        with sharedmem_pool(total_cores) as pool:
            def work(i):
                coeff = coeffs[i]
                thresh = threshold
                ne.evaluate('coeff * (real(abs(coeff)) >= thresh)',
                            out=thresh_coeffs[i])

            pool.map(work, range(self.redundancy))

        thresh_coeffs = np.array(thresh_coeffs)
        return self.inverse_transform(thresh_coeffs,
                                      real=real,
                                      do_norm=do_norm,
                                      total_cores=total_cores)

    def inverse_transform(self,
                          coeffs,
                          real=False,
                          do_norm=True,
                          total_cores=None):
        if total_cores is None:
            total_cores = self.total_cores

        # result = sm.full((self.height, self.width), 0, dtype='complex128')
        result = sm.full((self.height, self.width), 0, dtype=DTYPES[1])

        with sharedmem_pool(total_cores) as pool:
            def work(i):
                # local_fft = make_fft(self.width, self.height)
                # coeff_f = my_fft_shift(self._fft(coeffs[i], local_fft))

                self.add_inverse_part(coeffs[i], i, result, do_norm, pool)

                # coeff_fourier = my_fft_shift(self._fft(coeffs[i]))
                # assert coeff_fourier is not None  # silence pyflakes
                # spec = self._spectrograms[i]
                # dual_w = self.dual_frame_weight

                # temp = np.zeros((self.height, self.width),
                #                 dtype=DTYPES[1])

                # if do_norm:
                #     norm = self._space_norms[i]
                #     ne.evaluate('spec * norm * coeff_fourier / dual_w',
                #                 out=temp)
                # else:
                #     ne.evaluate('spec * coeff_fourier / dual_w', out=temp)

                # with pool.critical:
                #     result[:] += temp

            pool.map(work, list(range(len(self._spectrograms))))

        result = np.array(result)
        if real:
            # return np.real(ifft2(my_ifft_shift(result)))
            return np.real(self._ifft(my_ifft_shift(result)))
        else:
            # return ifft2(my_ifft_shift(result))
            return self._ifft(my_ifft_shift(result))

import numpy as np
import numexpr as ne

from healpy.pixelfunc import xyf2pix, npix2nside, reorder, nside2npix
from healpy.pixelfunc import pix2vec, vec2ang, get_interp_weights

from healpy.rotator import euler_matrix_new

from AlphaTransformParallel import sharedmem_pool

import sharedmem as sm

# NOTE: Preliminary tests show that no speedup is gained
#       (or in fact, a slowdone is achieved) when parallelizing the following:
#       * _blend_masks
#       * reassemble_faces
#       * get_faces
#       Therefore, these functions should be used in their 'old' version!

def make_index_par(nside, num_cores=None):
    # original sequential code:
    # index = np.array([xyf2pix(nside, x, range(nside), 0, True)
    #                   for x in range(nside-1, -1, -1)])

    index = sm.empty((nside, nside), dtype='int64')
    y = np.arange(nside)

    with sharedmem_pool(num_cores, numexpr=False) as pool:
        # we want to divide range(nside) into N chunks of (almost) equal size
        real_core_num = pool.np
        N = 10 * real_core_num # this can be changed if desired!
        chunk_size = 1 + nside // N # NOTE: N * chunk_size >= nside

        def work(i):
            lower = (N-(i+1))*chunk_size
            upper = (N-i)*chunk_size
            if lower < nside:
                index[lower:upper] = [xyf2pix(nside, x, y, 0, True)
                                      for x in range(lower, min(nside,upper))]

        pool.map(work, range(N-1,-1,-1))

    return np.array(index[::-1])

def get_faces_par(im, nested=False, num_cores=None):
    npix = np.shape(im)[0]
    assert npix % 12 == 0

    nside = npix2nside(npix)

    if not nested:
        new_im = reorder(im, r2n=True)
    else:
        new_im = im

    # index = np.array([xyf2pix(nside, x, range(nside), 0, True)
    #                     for x in range(nside-1, -1, -1)])
    index = make_index_par(nside, num_cores)

    CubeFace = sm.empty((12, nside, nside))
    with sharedmem_pool(num_cores, numexpr=False) as pool:
        # for face in range(12):
        #     CubeFace[face] = np.resize(new_im[index + nside**2 * face],
        #                                (nside, nside))
        def work(i):
            CubeFace[i] = np.resize(new_im[index + nside**2 * i],
                                    (nside, nside))

        pool.map(work, range(12))

    return np.array(CubeFace)


def get_faces(Imag, nested=False, num_cores=None):
    npix = np.shape(Imag)[0]
    assert npix % 12 == 0

    nside = npix2nside(npix)
    CubeFace = np.zeros((12, nside, nside))

    if not nested:
        NewIm = reorder(Imag, r2n=True)
    else:
        NewIm = Imag

    # index = np.array([xyf2pix(nside, x, range(nside), 0, True)
    #                     for x in range(nside-1, -1, -1)
    #                 ])
    index = make_index_par(nside, num_cores)

    for face in range(12):
        CubeFace[face] = np.resize(NewIm[index + nside**2 * face], (nside, nside))

    return CubeFace


def set_faces_par(CubeFace, nested=False, num_cores=None):
    npix = np.size(CubeFace)
    assert npix % 12 == 0

    nside = npix2nside(npix)

    index = make_index_par(nside, num_cores)
    # index = np.array([xyf2pix(nside, x, range(nside), 0, True)
    #                     for x in range(nside-1, -1, -1)])

    imag = sm.empty((npix))
    with sharedmem_pool(num_cores, numexpr=False) as pool:
        def work(i):
             for face in range(12):
                 imag[index + nside**2 * face] = np.resize(CubeFace[face],
                                                           (nside, nside))
            #imag[index + nside**2 * face] = np.resize(CubeFace[face],
            #                                          (nside,nside))

        pool.map(work, range(12))

    imag = np.array(imag)

    if not nested:
        new_im = reorder(imag, n2r=True)
    else:
        new_im = imag

    return new_im

def set_faces(CubeFace, nested=False):
    npix = np.size(CubeFace)
    assert npix % 12 == 0

    nside = npix2nside(npix)
    Imag = np.zeros((npix))

    index = np.array([xyf2pix(nside, x, range(nside), 0, True)
                        for x in range(nside-1, -1, -1)
                    ])

    for face in range(12):
        Imag[index + nside**2 * face] = np.resize(CubeFace[face], (nside, nside))

    if not nested:
        NewIm = reorder(Imag, n2r=True)
    else:
        NewIm = Imag

    return NewIm

def rotate_map(Imag,a1,a2,a3,X=True,Y=False,ZYX=False,deg=False,nested=False):
    # X :   rotation a1 around original Z
    #       rotation a2 around interm   X
    #       rotation a3 around final    Z
    #            DEFAULT,  classical mechanics convention

    #  Y :  rotation a1 around original Z
    #       rotation a2 around interm   Y
    #       rotation a3 around final    Z
    #            quantum mechanics convention (override X)

    #  ZYX :rotation a1 around original Z
    #       rotation a2 around interm   Y
    #       rotation a3 around final    X
    #            aeronautics convention (override X)
    #  * these last three keywords are obviously mutually exclusive *

    npix=np.shape(Imag)[0]
    nside=npix2nside(npix)
    indices=np.arange(0,npix)
    ang_coord=pix2vec(nside, indices,nested)
    ang_coord_array=np.vstack((ang_coord[0],ang_coord[1],ang_coord[2]))
    eul=euler_matrix_new(a1,a2,a3,X=X,Y=Y,ZYX=ZYX,deg=deg)
    new_coord=np.dot(eul,ang_coord_array)
    theta_arr,phi_arr=vec2ang(new_coord.T)
    neigh,weigh=get_interp_weights(nside,theta_arr,phi=phi_arr,nest=nested,lonlat=False)
    thr_val=1e-8
    weigh[np.where(np.abs(weigh)<thr_val)]=0
    weigh=weigh/np.sum(weigh,axis=0)
    rotIm=np.zeros_like(Imag)
    for k in range(neigh.shape[0]):
        rotIm=rotIm+weigh[k]*Imag[neigh[k]]
    return rotIm

def rotate_map_inv(Imag,a1,a2,a3,X=True,Y=False,ZYX=False,deg=False,nested=False):
    npix=np.shape(Imag)[0]
    nside=npix2nside(npix)
    indices=np.arange(0,npix)
    ang_coord=pix2vec(nside, indices,nested)
    ang_coord_array=np.vstack((ang_coord[0],ang_coord[1],ang_coord[2]))
    eul=np.linalg.inv(euler_matrix_new(a1,a2,a3,X=X,Y=Y,ZYX=ZYX,deg=deg))
    new_coord=np.dot(eul,ang_coord_array)
    theta_arr,phi_arr=vec2ang(new_coord.T)
    neigh,weigh=get_interp_weights(nside,theta_arr,phi=phi_arr,nest=nested,lonlat=False)
    thr_val=1e-8
    weigh[np.where(np.abs(weigh)<thr_val)]=0
    weigh=weigh/np.sum(weigh,axis=0)
    rotIm=np.zeros_like(Imag)
    for k in range(neigh.shape[0]):
        rotIm=rotIm+weigh[k]*Imag[neigh[k]]
    return rotIm

def compose_polar_faces(faces, overlap, num_cores=None):
    nside = faces[0].shape[0]

    assert faces[0].shape[1] == nside
    assert overlap <= nside//4

    # FIXME: If this is changed, change the constructor of PW_Transform!
    nside_output = 2 * (3*nside//4 + overlap)

    cut = nside//4 - overlap
    size = nside - cut

    # north = np.ones((4, size, size))
    north = sm.empty((4, size, size))

    # north[0] = np.rot90(faces[1][:size, cut:], -1)
    # north[1] = np.rot90(faces[2][:size, cut:], -2)
    # north[2] = np.rot90(faces[0][:size, cut:], 0)
    # north[3] = np.rot90(faces[3][:size, cut:], -3)

    first_index = (1,2,0,3)
    second_index = (-1,-2,0,-3)

    with sharedmem_pool(num_cores, numexpr=False) as pool:
        def work(i):
            north[i] = np.rot90(faces[first_index[i]][:size, cut:],
                                second_index[i])

        pool.map(work, range(4))

    output = np.ones((2, nside_output, nside_output))

    output[0] = np.vstack(( np.hstack((north[0], north[1])),
                            np.hstack((north[2], north[3]))
                            ))
    del north

    # south = np.ones((4, size, size))
    south = sm.empty((4, size, size))

    first_index = (9,8,10,11)
    second_index = (1,0,2,3)

    # south[0] = np.rot90(faces[9][cut:, :size], 1)
    # south[1] = np.rot90(faces[8][cut:, :size], 0)
    # south[2] = np.rot90(faces[10][cut:, :size], 2)
    # south[3] = np.rot90(faces[11][cut:, :size], 3)

    with sharedmem_pool(num_cores, numexpr=False) as pool:
        def work(i):
            south[i] = np.rot90(faces[first_index[i]][cut:, :size],
                                second_index[i])

        pool.map(work, range(4))

    output[1] = np.vstack(( np.hstack((south[0], south[1])),
                            np.hstack((south[2], south[3]))
                            ))

    return output

def compose_polar_faces_old(faces, overlap):
    nside = faces[0].shape[0]

    assert faces[0].shape[1] == nside

    assert overlap <= nside//4

    # FIXME: If this is changed, change the constructor of PW_Transform!
    nside_output = 2 * (3*nside//4 + overlap)

    cut = nside//4 - overlap
    size = nside - cut

    north = np.ones((4, size, size))

    north[0] = np.rot90(faces[1][:size, cut:], -1)
    north[1] = np.rot90(faces[2][:size, cut:], -2)
    north[2] = np.rot90(faces[0][:size, cut:], 0)
    north[3] = np.rot90(faces[3][:size, cut:], -3)

    output = np.ones((2, nside_output, nside_output))

    output[0] = np.vstack(( np.hstack((north[0], north[1])),
                            np.hstack((north[2], north[3]))
                            ))

    south = np.ones((4, size, size))

    south[0] = np.rot90(faces[9][cut:, :size], 1)
    south[1] = np.rot90(faces[8][cut:, :size], 0)
    south[2] = np.rot90(faces[10][cut:, :size], 2)
    south[3] = np.rot90(faces[11][cut:, :size], 3)


    output[1] = np.vstack(( np.hstack((south[0], south[1])),
                            np.hstack((south[2], south[3]))
                            ))

    return output


def _face_masks_old(size):

    def _cos_step(x, size):
        return ((0 <= x) * (-0.5 * np.cos(np.pi*x/size) + 0.5) * (x < size)
                + (x >= size))

    masks = np.zeros((4, size, size))

    # one of the values is shifted to avoid for the sum of two masks to vanish
    x, y = np.meshgrid(np.linspace(0, size-1, size), np.linspace(1, size, size))

    masks[0] = np.rot90(_cos_step(size-x, size) * _cos_step(2*y-x, size), 1)
    masks[1] = np.rot90(_cos_step(size-y+1, size) * _cos_step(2*x-y+2, size), 1)
    masks[2] = np.rot90(masks[0], 2)
    masks[3] = np.rot90(masks[1], 2)

    return masks

def _face_masks(size, num_cores=None):

    def _cos_step(x, size):
        pi = np.pi
        return ne.evaluate('(0 <= x) * (-0.5 * cos(pi*x/size) + 0.5)'
                           '* (x < size) + (x >= size)')
        # return ((0 <= x) * (-0.5 * np.cos(np.pi*x/size) + 0.5) * (x < size)
        #         + (x >= size))

    # masks = np.zeros((4, size, size))
    masks = sm.empty((4, size, size))

    # one of the values is shifted to avoid for the sum of two masks to vanish
    x, y = np.meshgrid(np.linspace(0, size-1, size),
                       np.linspace(1, size, size))

    # masks[0] = np.rot90(_cos_step(size-x, size) * _cos_step(2*y-x, size), 1)
    # masks[1] = np.rot90(_cos_step(size-y+1, size) * _cos_step(2*x-y+2, size), 1)
    # masks[2] = np.rot90(masks[0], 2)
    # masks[3] = np.rot90(masks[1], 2)

    with sharedmem_pool(num_cores, numexpr=True) as pool:
        def work(i):
            if i == 0:
                masks[0] = np.rot90(_cos_step(size-x, size) * _cos_step(2*y-x, size), 1)
                masks[2] = np.rot90(masks[0], 2)
            else:
                masks[1] = np.rot90(_cos_step(size-y+1, size) * _cos_step(2*x-y+2, size), 1)
                masks[3] = np.rot90(masks[1], 2)

        pool.map(work, range(2))

    return np.array(masks)


def compose_equatorial_faces_old(face, overlap):
    nside = face[0].shape[0]

    assert face[0].shape[1] == nside

    assert overlap <= 3*nside//4

    extra = nside//4 + overlap

    nside_output = nside + 2*extra

    output = np.ones((4, nside_output, nside_output))

    masks = _face_masks_old(extra)

    for i in range(4):
        rows = []
        temp_row = []

        temp_row.append(face[4+(i+1)%4][-extra:, -extra:])
        temp_row.append(face[i][-extra:, :])
        temp_row.append(( masks[0] * np.rot90(face[i][-extra:, -extra:], -1)
                        + masks[1] * np.rot90(face[(i+3)%4][:extra, :extra], 1)
                        ) / (masks[0] + masks[1]))

        rows.append(np.hstack(temp_row))

        temp_row = []

        temp_row.append(face[8+i][:, -extra:])
        temp_row.append(face[4+i])
        temp_row.append(face[(i+3)%4][:, :extra])

        rows.append(np.hstack(temp_row))

        temp_row = []

        temp_row.append(( masks[2] * np.rot90(face[8+(i+3)%4][:extra, :extra], -1)
                        + masks[3] * np.rot90(face[8+i][-extra:, -extra:], 1)
                        ) / (masks[2] + masks[3]))

        temp_row.append(face[8+(i+3)%4][:extra, :])
        temp_row.append(face[4+(i+3)%4][:extra, :extra])

        rows.append(np.hstack(temp_row))

        output[i] = np.vstack(rows)

    return output

def compose_equatorial_faces(faces, overlap, num_cores=None):
    nside = faces[0].shape[0]

    assert faces[0].shape[1] == nside
    assert overlap <= 3*nside//4

    extra = nside//4 + overlap
    nside_output = nside + 2*extra
    masks = _face_masks(extra, num_cores)

    output = sm.empty((4, nside_output, nside_output))

    with sharedmem_pool(num_cores, numexpr=False) as pool:
        def work(i):
            rows = []
            temp = []

            temp.append(faces[4+(i+1)%4][-extra:, -extra:])
            temp.append(faces[i][-extra:, :])
            temp.append((masks[0] * np.rot90(faces[i][-extra:, -extra:], -1)
                         + masks[1]*np.rot90(faces[(i+3)%4][:extra,:extra], 1))
                        / (masks[0] + masks[1]))

            rows.append(np.hstack(temp))

            temp = []

            temp.append(faces[8+i][:, -extra:])
            temp.append(faces[4+i])
            temp.append(faces[(i+3)%4][:, :extra])

            rows.append(np.hstack(temp))

            temp = []

            temp.append((masks[2]*np.rot90(faces[8+(i+3)%4][:extra,:extra], -1)
                         + masks[3] * np.rot90(faces[8+i][-extra:, -extra:], 1))
                        / (masks[2] + masks[3]))

            temp.append(faces[8+(i+3)%4][:extra, :])
            temp.append(faces[4+(i+3)%4][:extra, :extra])

            rows.append(np.hstack(temp))

            output[i] = np.vstack(rows)

        pool.map(work, range(4))

    return np.array(output)

def _blend_masks(nside, size):

    def _cos(x):
        return -0.5 * np.cos(np.pi*x) + 0.5

    x, y = np.meshgrid( np.linspace(1, nside, nside),
                        np.hstack([ np.zeros(3*nside//4-size),
                                    np.linspace(0, 1, size),
                                    np.ones(nside//4)
                                    ])
                        )

    blend_masks = np.zeros((8, nside, nside))

    blend_masks[0] = _cos(y)
    blend_masks[1] = np.rot90(blend_masks[0], -1)
    blend_masks[2] = np.rot90(blend_masks[0], -2)
    blend_masks[3] = np.rot90(blend_masks[0], -3)

    blend_masks[4] = blend_masks[0] * blend_masks[3]
    blend_masks[5] = blend_masks[1] * blend_masks[2]

    blend_masks[6] = (1 - blend_masks[0]) * (1 - blend_masks[1])
    blend_masks[7] = (1 - blend_masks[2]) * (1 - blend_masks[3])

    return blend_masks


# transition in [north, east, south, west, north-west, south-east,
#               north-polar, south-polar] direction at polar regions
def _blend_masks_par(nside, size, num_cores=None):

    def _cos(x):
        # return -0.5 * np.cos(np.pi*x) + 0.5
        pi = np.pi
        return ne.evaluate('-0.5 * cos(pi * x) + 0.5')

    x, y = np.meshgrid(np.linspace(1, nside, nside),
                       np.hstack([np.zeros(3*nside//4-size),
                                  np.linspace(0, 1, size),
                                  np.ones(nside//4)]))

    blend_masks = sm.empty((8, nside, nside))

    blend_masks[0] = _cos(y)
    # blend_masks[1] = np.rot90(blend_masks[0], -1)
    # blend_masks[2] = np.rot90(blend_masks[0], -2)
    # blend_masks[3] = np.rot90(blend_masks[0], -3)
    with sharedmem_pool(num_cores, numexpr=False) as pool:
        def work(i):
            blend_masks[i] = np.rot90(blend_masks[0], -i)

        pool.map(work, range(1,4))

    # blend_masks[4] = blend_masks[0] * blend_masks[3]
    # blend_masks[5] = blend_masks[1] * blend_masks[2]
    # blend_masks[6] = (1 - blend_masks[0]) * (1 - blend_masks[1])
    # blend_masks[7] = (1 - blend_masks[2]) * (1 - blend_masks[3])
    first_index = (0, 1, 0, 2)
    second_index = (3, 2, 1, 3)
    with sharedmem_pool(num_cores) as pool:
        def work(i):
            if i in {4,5}:
                blend_masks[i] = (blend_masks[first_index[i-4]] *
                                  blend_masks[second_index[i-4]])
            elif i in {6,7}:
                blend_masks[i] = ((1 - blend_masks[first_index[i-4]]) *
                                  (1 - blend_masks[second_index[i-4]]))

        pool.map(work, range(4,8))

    return np.array(blend_masks)

def reassemble_faces_par(nside,
                         polar_faces,
                         eq_faces,
                         margin,
                         transition,
                         num_cores=None):
    out = sm.empty((12, nside, nside))

    blend_masks = _blend_masks_par(nside, transition, num_cores)

    offset = 3*nside//4 + transition + margin
    size = offset - margin

    with sharedmem_pool(num_cores, numexpr=False) as pool:
        def work(i):
            if i <= 3:
                ### north polar faces
                out[i][:size, -size:] = np.rot90(polar_faces[0],
                                                 i)[offset:-margin,
                                                    margin:offset]
                out[i] *= blend_masks[6]
            else:
                ### south polar faces
                j = i - 8
                out[i][-size:, :size] = np.rot90(polar_faces[1],
                                                 -j)[margin:offset,
                                                     offset:-margin]
                out[i] *= blend_masks[7]

        pool.map(work, [0,1,2,3,8,9,10,11])

    extra = nside//4 + transition
    offset = extra + margin

    with sharedmem_pool(num_cores, numexpr=False) as pool:
        def work(i):
            if i < 4:
                # north
                out[i][-extra:] += (blend_masks[0][-extra:]
                                    * eq_faces[i][margin:offset,
                                                  offset:-offset])
                # east
                out[i][:, :extra] += (blend_masks[1][:, :extra]
                                      * eq_faces[(i+1)%4][offset:-offset,
                                                          -offset:-margin])
            elif i < 8:
                # north-west
                out[i][-extra:, -extra:] += (blend_masks[4][-extra:, -extra:]
                                             * eq_faces[(i-1)%4][margin:offset,
                                                                 margin:offset])
                # middle
                out[i] += (eq_faces[i-4][offset:-offset, offset:-offset])
                # south-east
                out[i][:extra, :extra] += (blend_masks[5][:extra, :extra]
                                           * eq_faces[(i-3)%4][-offset:-margin,
                                                               -offset:-margin])
            else:
                # west
                out[i][:, -extra:] += (blend_masks[3][:, -extra:]
                                       * eq_faces[i-8][offset:-offset,
                                                       margin:offset])
                # south
                out[i][:extra] += (blend_masks[2][:extra]
                                   * eq_faces[(i-7)%4][-offset:-margin,
                                                       offset:-offset])

        pool.map(work, range(12))

    with sharedmem_pool(num_cores, numexpr=False) as pool:
        def work(i):
            if i < 4:
                out[i] /= (blend_masks[6] + blend_masks[0] + blend_masks[1])
            elif i < 8:
                out[i] /= (1 + blend_masks[4] + blend_masks[5])
            else:
                out[i] /= (blend_masks[7] + blend_masks[2] + blend_masks[3])

        pool.map(work, range(12))

    return np.array(out)

def reassemble_faces(nside, polar_faces, equatorial_faces, margin, transition):
    output = np.zeros((12, nside, nside))

    blend_masks = _blend_masks(nside, transition)

    offset = 3*nside//4 + transition + margin

    size = offset - margin

    for i in range(4):
        ### north polar faces
        output[i][:size, -size:] = np.rot90(polar_faces[0], i)[offset:-margin,
                                                               margin:offset]
        output[i] *= blend_masks[6]
        ### south polar faces
        output[i+8][-size:, :size] = np.rot90(polar_faces[1], -i)[margin:offset,
                                                                  offset:-margin]
        output[i+8] *= blend_masks[7]

    ### equatorial faces and overlapping areas

    extra = nside//4 + transition

    offset = extra + margin

    for i in range(4):
        # north-west
        output[4+(i+1)%4][-extra:, -extra:] += (blend_masks[4][-extra:, -extra:]
                        * equatorial_faces[i][margin:offset, margin:offset])
        # middle
        output[4+i] += (equatorial_faces[i][offset:-offset, offset:-offset])
        # south-east
        output[4+(i+3)%4][:extra, :extra] += (blend_masks[5][:extra, :extra]
                        * equatorial_faces[i][-offset:-margin, -offset:-margin])

        # west
        output[8+i][:, -extra:] += (blend_masks[3][:, -extra:]
                        * equatorial_faces[i][offset:-offset, margin:offset])
        # south
        output[8+(i+3)%4][:extra] += (blend_masks[2][:extra]
                        * equatorial_faces[i][-offset:-margin, offset:-offset])

        # north
        output[i][-extra:] += (blend_masks[0][-extra:]
                        * equatorial_faces[i][margin:offset, offset:-offset])
        # east
        output[(i+3)%4][:, :extra] += (blend_masks[1][:, :extra]
                        * equatorial_faces[i][offset:-offset, -offset:-margin])

    for i in range(4):
        output[i] /= (blend_masks[6] + blend_masks[0] + blend_masks[1])
        output[i+8] /= (blend_masks[7] + blend_masks[2] + blend_masks[3])
        output[i+4] /= (1 + blend_masks[4] + blend_masks[5])

    return output

# NOTE: I did not try to parallelize this function since it is just too horrible
#       Also, I do not expect the performance gain to be huge!
def arrange_faces(faces, overlap=0, v=True, h=True, orientation=False):
    nside = faces[0].shape[0]
    assert nside == faces[0].shape[1]
    #assert nside2npix(nside) == np.size(faces)

    if orientation:
        faces = faces.copy()
        for i in range(12):
            faces[i] = np.rot90(faces[i], -1)

    temp = -float('inf') * np.ones((32, nside, nside))

    temp[1] = faces[11]
    temp[2] = faces[7]
    temp[3] = faces[2]
    temp[10] = faces[10]
    temp[11] = faces[6]
    temp[12] = faces[1]
    temp[19] = faces[9]
    temp[20] = faces[5]
    temp[21] = faces[0]
    temp[28] = faces[8]
    temp[29] = faces[4]
    temp[30] = faces[3]

    if not (overlap > 0):
        return np.vstack([np.hstack(temp[0:8]),
                          np.hstack(temp[8:16]),
                          np.hstack(temp[16:24]),
                          np.hstack(temp[24:36])])

    assert overlap % 1 == 0
    masks = _face_masks(overlap)

    start = nside - overlap

    if h:
        if v:
            temp[0][:overlap, start:] = (
                  masks[2] * np.rot90(faces[11][:overlap, :overlap], -1)
                + masks[3] * np.rot90(faces[8][start:, start:], 1)
                                        ) / (masks[2] + masks[3])
            temp[9][:overlap, start:] = (
                  masks[2] * np.rot90(faces[10][:overlap, :overlap], -1)
                + masks[3] * np.rot90(faces[11][start:, start:], 1)
                                        ) / (masks[2] + masks[3])
            temp[18][:overlap, start:] = (
                  masks[2] * np.rot90(faces[9][:overlap, :overlap], -1)
                + masks[3] * np.rot90(faces[10][start:, start:], 1)
                                        ) / (masks[2] + masks[3])
            temp[27][:overlap, start:] = (
                  masks[2] * np.rot90(faces[8][:overlap, :overlap], -1)
                + masks[3] * np.rot90(faces[9][start:, start:], 1)
                                        )/(masks[2] + masks[3])

            temp[4][start:, :overlap] = (
                  masks[0] * np.rot90(faces[2][start:, start:], -1)
                + masks[1] * np.rot90(faces[1][:overlap, :overlap], 1)
                                            ) / (masks[0] + masks[1])
            temp[13][start:, :overlap] = (
                  masks[0] * np.rot90(faces[1][start:, start:], -1)
                + masks[1] * np.rot90(faces[0][:overlap, :overlap], 1)
                                            ) / (masks[0] + masks[1])
            temp[22][start:, :overlap] = (
                  masks[0] * np.rot90(faces[0][start:, start:], -1)
                + masks[1] * np.rot90(faces[3][:overlap, :overlap], 1)
                                            ) / (masks[0] + masks[1])
            temp[31][start:, :overlap] = (
                  masks[0] * np.rot90(faces[3][start:, start:], -1)
                + masks[1] * np.rot90(faces[2][:overlap, :overlap], 1)
                                            ) / (masks[0] + masks[1])
        else:
            temp[0][:overlap, start:] = (
                    masks[3] * np.rot90(faces[8][start:, start:], 1))
            temp[9][:overlap, start:] = (
                    masks[3] * np.rot90(faces[11][start:, start:], 1))
            temp[18][:overlap, start:] = (
                    masks[3] * np.rot90(faces[10][start:, start:], 1))
            temp[27][:overlap, start:] = (
                    masks[3] * np.rot90(faces[9][start:, start:], 1))

            temp[4][start:, :overlap] = (
                    masks[1] * np.rot90(faces[1][:overlap, :overlap], 1))
            temp[13][start:, :overlap] = (
                    masks[1] * np.rot90(faces[0][:overlap, :overlap], 1))
            temp[22][start:, :overlap] = (
                    masks[1] * np.rot90(faces[3][:overlap, :overlap], 1))
            temp[31][start:, :overlap] = (
                    masks[1] * np.rot90(faces[2][:overlap, :overlap], 1))
    elif v:
        temp[0][:overlap, start:] = (
                    masks[2] * np.rot90(faces[11][:overlap, :overlap], -1))
        temp[9][:overlap, start:] = (
                    masks[2] * np.rot90(faces[10][:overlap, :overlap], -1))
        temp[18][:overlap, start:] = (
                    masks[2] * np.rot90(faces[9][:overlap, :overlap], -1))
        temp[27][:overlap, start:] = (
                    masks[2] * np.rot90(faces[8][:overlap, :overlap], -1))

        temp[4][start:, :overlap] = (
                    masks[0] * np.rot90(faces[2][start:, start:], -1))
        temp[13][start:, :overlap] = (
                    masks[0] * np.rot90(faces[1][start:, start:], -1))
        temp[22][start:, :overlap] = (
                    masks[0] * np.rot90(faces[0][start:, start:], -1))
        temp[31][start:, :overlap] = (
                    masks[0] * np.rot90(faces[3][start:, start:], -1))

    return np.vstack([np.hstack(temp[0:8]),
                      np.hstack(temp[8:16]),
                      np.hstack(temp[16:24]),
                      np.hstack(temp[24:36])])

def normalize(image):
    minval = np.min(image)
    output = image - minval
    maxval = np.max(output)
    return [output/maxval, minval, maxval]

def denormalize(image, minval, maxval):
    return image*maxval + minval

def l2_normalize(image):
    minv = np.min(image)
    output = image - minv
    norm = np.linalg.norm(output)
    return [output/norm, minv, norm]

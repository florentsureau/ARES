import math
import healpy as hp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from healpy.rotator import *
from astropy.io import fits as fits
import copy
from os.path import expanduser
import sys

DtoR=math.pi/180.

def getAllIndicesForFace(nside,face,nested=False):
    #Note X: column, y : Rows, Row-Major ordering: Face[y,x]
    index=np.array([hp.xyf2pix(nside,range(nside),y,face,nest=nested)
                                                        for y in range(nside)])
    return np.reshape(index,(nside,nside))

def getAllFacesIndex(nside,nested=False):
    #RETURN FOR EACH PIXEL ITS INDEX IN CUBE OF FACES
    if not hp.isnsideok(nside):
        raise ValueError('Incorrect nside')
    npix=hp.nside2npix(nside)
    index_map=np.zeros((npix))
    offset=nside*nside
    #Number one face
    if (nested):
        #Only need to do it for one face
        index=np.array([hp.xyf2pix(nside,range(nside),y,0,nested)
                                                        for y in range(nside)])
        for face in range(12):
            index_map[index+face*offset]=np.linspace(0,offset-1,offset)+face*offset
    else:
        for face in range(12):
            index=np.array([hp.xyf2pix(nside,range(nside),y,face,nested)
                                                        for y in range(nside)])
            index_map[index]=np.linspace(0,offset-1,offset)+face*offset
    return index_map

def get_face(Imag,face,nested=False):
    #Get a given face in the image
    npix=np.shape(Imag)[0]
    if not hp.isnpixok(npix):
        raise ValueError('Incorrect map size {0}'.format(npix))
    nside = hp.npix2nside(npix)
    taille_face = npix/12
    cote=int(math.sqrt(taille_face))
    CubeFace = np.zeros((cote,cote))
    if (nested!=True):
        NewIm=hp.reorder(Imag, r2n = True)
    else :
        NewIm=Imag
    index = np.zeros((cote,cote))
    index=getAllIndicesForFace(nside,0,nested=True)
    print("Process Face {0}".format(face))
    CubeFace[:,:] =np.resize(NewIm[index+taille_face*face],(cote,cote))
    return CubeFace


def get_all_faces(Imag,nested=False):
    #Get the cube of faces
    npix=np.shape(Imag)[0]
    if not hp.isnpixok(npix):
        raise ValueError('Incorrect map size {0}'.format(npix))
    nside =  hp.npix2nside(npix)
    taille_face = npix/12
    cote=int(math.sqrt(taille_face))
    print(cote)
    CubeFace = np.zeros((12,cote,cote))
    if (nested!=True):
        NewIm=hp.reorder(Imag, r2n = True)
    else :
        NewIm=Imag
    index=np.zeros((cote,cote))
    index=getAllIndicesForFace(nside,0,nested=True)
    for face in range(12):
        print("Process Face {0}".format(face))
        CubeFace[face,:,:] =np.resize(NewIm[index+taille_face*face],(cote,cote))
        #plt.figure(),imshow(np.log10(1+CubeFace[face,:,:]*1e6))
        #plt.title("face {0}".format(face)),plt.colorbar()
    return CubeFace

def put_all_faces(CubeFace,nested=False):
    #From a cube of face, reconstruct an HEALPix image
    npix=np.size(CubeFace)
    if not hp.isnpixok(npix):
        raise ValueError('Incorrect cube size {0}'.format(npix))
    nside = hp.npix2nside(npix)
    taille_face = npix/12
    cote=int(math.sqrt(taille_face))
    Imag = np.zeros((npix))
    index = np.zeros((cote,cote))
    index=getAllIndicesForFace(nside,0,nested=True)
    for face in range(12):
        print("Process Face {0}".format(face))
        Imag[index+taille_face*face]=np.resize(CubeFace[face,:,:],(cote,cote))

    if (nested!=True):
        NewIm=hp.reorder(Imag, n2r = True)
    else:
        NewIm=Imag
    return NewIm

def constructExtendedFaceIndices(nside,face,PatchWidth,nested=False,Verbose=False):
    #THE CATCH IS THAT THE NEIGHBOUR ORDERING IS NOT WITH RESPECT TO THE CURRENT
    #PIXEL, BUT TO A POINT OUTSIDE THE SPHERE SO FOR A PIXEL IN A CAP WE HAVE
    # TO IDENTIFY:
    # 1) SINGULAR PIXELS, WITH ONLY 7-NEIGHBOUR
    # 2) WHAT IS THE RELEVANT ORIENTATION WITH RESPECT TO THE CURRENT PIXEL
    #    IN THE GRID
    #(1) IS DONE BY LOOKING AT A VALUE OF -1 IN THE LIST OF NEIGHBOUR PIXELS
    #(2) IS DONE BY TAKING A REFERENCE POINT IN THE GRID WITH RESPECT TO CURRENT
    #    PIXEL (UPPER PIXEL FOR RIGHT BORDER EXTENSION, AND RIGHT PIXEL FOR TOP
    #    BORDER EXTENSION) AND FILLING THE GRID ACCORDINGLY
    FaceExtension=np.zeros((nside+PatchWidth-1,nside+PatchWidth-1))-1
    FaceExtension[PatchWidth-1:nside+PatchWidth-1,0:nside]=getAllIndicesForFace(
                                    nside,face,nested=nested).astype('int64')

    #RIGHT BORDER EXTENSION
    if not hp.isnsideok(nside):
        raise ValueError('Incorrect map nside {0}'.format(nside))
    for kx in range(nside,nside+PatchWidth-1):
        lstPix=FaceExtension[PatchWidth-1:nside+PatchWidth-1,kx-1].astype(
                                                              'int64').flatten()
        PixBorder=np.array(hp.get_all_neighbours(nside,lstPix,nest=nested))
        for ky in range(PatchWidth-1,nside+PatchWidth-2):
            UpPix=FaceExtension[ky+1,kx-1]
            if(-1 in PixBorder[0:8,ky-PatchWidth+1]):
                if(Verbose):
                    print("({0},{1}):".format(kx,ky))
                    print("CANNOT PROCESS SINGULAR 7-NEIGHBOUR PIXEL")
            else:
                startInd=(np.where(PixBorder[0:8,ky-PatchWidth+1]==UpPix))[0][0]
                FaceExtension[ky+1,kx]=PixBorder[(startInd+1)%8,ky-PatchWidth+1]
                FaceExtension[ky,kx]=PixBorder[(startInd+2)%8,ky-PatchWidth+1]
                if(ky>0) and (FaceExtension[ky-1,kx]==-1):
                    FaceExtension[ky-1,kx]=PixBorder[(startInd+3)%8,
                                                                ky-PatchWidth+1]

    #BOTTOM BORDER EXTENSION
    for ky in range(PatchWidth-2,-1,-1):
        lstPix=FaceExtension[ky+1,0:nside].astype('int64').flatten()
        PixBorder=np.array(hp.get_all_neighbours(nside,lstPix,nest=False))
        for kx in range(1,nside):
            LeftPix=FaceExtension[ky+1,kx-1]
            if(-1 in PixBorder[0:8,kx]):
                if(Verbose):
                    print("({0},{1}):".format(kx,ky))
                    print("CANNOT PROCESS SINGULAR 7-NEIGHBOUR PIXEL")
            else:
                lstStart=(np.where(PixBorder[0:8,kx]==LeftPix))
                if(len(lstStart[0])>0):
                    startInd=lstStart[0][0]
                    if(FaceExtension[ky,kx-1]==-1):
                        FaceExtension[ky,kx-1]=PixBorder[(startInd-1)%8,kx]
                    FaceExtension[ky,kx]=PixBorder[(startInd-2)%8,kx]
                    FaceExtension[ky,kx+1]=PixBorder[(startInd-3)%8,kx]

    #DIAGONAL EXTENSION FROM RIGHT BORDER EXTENSION
    for ky in range(PatchWidth-2,-1,-1):
        lstPix=FaceExtension[ky+1,nside:nside+PatchWidth-1].astype(
                                                            'int64').flatten()
        if(-1 in lstPix):
            if(Verbose):
                print("NOTHING TO DO")
        else:
            PixBorder=np.array(hp.get_all_neighbours(nside,lstPix,nest=nested))
            for kx in range(nside,nside+PatchWidth-1):
                UpperPix=FaceExtension[ky+2,kx]
                if(-1 in PixBorder[0:8,kx-nside]):
                    if(Verbose):
                        print("CANNOT PROCESS PIXEL ({0},{1})".format(kx,ky))
                else:
                    lstStart=(np.where(PixBorder[0:8,kx-nside]==UpperPix))
                    if(len(lstStart[0])>0):
                        if(Verbose):
                            print("PROCESS PIXEL ({0},{1})".format(kx,ky))
                        startInd=lstStart[0][0]
                        FaceExtension[ky,kx]=PixBorder[(startInd+4)%8,kx-nside]
                        if(kx<nside+PatchWidth-2):
                            FaceExtension[ky,kx+1]=PixBorder[(startInd+5)%8,
                                                                    kx-nside]
    if Verbose:
        print(FaceExtension.astype('int64'))
    return FaceExtension.astype('int64')

def getAll2DPatches(Imag,PatchWidth,nested=False,BorderCap=False,Verbose=False):
    #From one image, get all patches
    npix=np.shape(Imag)[0]
    if not hp.isnpixok(npix):
        raise ValueError('Incorrect map size {0}'.format(npix))
    nside=hp.npix2nside(npix)
    if(nside<PatchWidth):
        raise ValueError('Incorrect PatchWidth {0} vs nside {1}'.format(
                                                            PatchWidth,nside))
    Patches2D=np.zeros((npix,PatchWidth,PatchWidth))
    NPatchPerFace=np.zeros((12),dtype='int64')
    NPatchPerSide=nside
    offsetPatchIndex=0
    #Process all patches without -1
    for f in range(12):
        if(Verbose):
            print("PROCESS FACE {0}".format(f))
        FaceIndices=constructExtendedFaceIndices(nside,f,PatchWidth,
                                                                nested=nested)
        ExtendedFace=Imag[FaceIndices.astype('int64')]
        FacePatches2D=np.array([ExtendedFace[y:y+PatchWidth,x:x+PatchWidth]
            for y in range(nside) for x in range(nside)
                    if -1 not in FaceIndices[y:y+PatchWidth,x:x+PatchWidth]])
        NPatchPerFace[f]=np.shape(FacePatches2D)[0]
        Patches2D[offsetPatchIndex:offsetPatchIndex+NPatchPerFace[f]]=\
                                                                FacePatches2D
        offsetPatchIndex=offsetPatchIndex+NPatchPerFace[f]
    Patches2D=Patches2D[0:offsetPatchIndex,:,:]
    return np.array(Patches2D),NPatchPerFace


def get2DPatchesPositions(nside,PatchWidth,PatchIndices,nested=False,
                                                                Verbose=False):
    #Get all pixel position on the sphere for the different patches
    if not hp.isnsideok(nside):
        raise ValueError('Incorrect map nside {0}'.format(nside))
    nside2=nside*nside
    NpatchesCap=(nside2-(PatchWidth-1)**2) #For Northern/Southern cap Patches
    NpatchesEq=nside2 # Number of Equatorial Patches
    npatches=NpatchesCap*8+NpatchesEq*4

    #Start by sorting the indices to process face by face
    SortingIndices=np.argsort(PatchIndices)
    sortedIndices=PatchIndices[SortingIndices]
    AllPatchesPosition=np.zeros((npatches,PatchWidth*PatchWidth),dtype='int64')

    NPatchPerFace=np.zeros((12),dtype='int64')
    NPatchPerSide=nside
    offsetPatchIndex=0
    #Process all patches without -1
    for f in range(12):
        if(Verbose):
            print("PROCESS FACE {0}".format(f))
        FaceIndices=constructExtendedFaceIndices(nside,f,PatchWidth,
                                                            nested=nested)
        FacePatches2D=np.array([np.ndarray.flatten(FaceIndices[
                    y:y+PatchWidth,x:x+PatchWidth])
                    for y in range(nside) for x in range(nside)
                    if -1 not in FaceIndices[y:y+PatchWidth,x:x+PatchWidth]])
        NPatchPerFace[f]=np.shape(FacePatches2D)[0]
        AllPatchesPosition[offsetPatchIndex:offsetPatchIndex+NPatchPerFace[f]]=\
                                                                FacePatches2D
        offsetPatchIndex=offsetPatchIndex+NPatchPerFace[f]
    PatchesPosition=AllPatchesPosition[PatchIndices]
    return PatchesPosition

def get2DPatchesCoverage(nside,PatchWidth,PatchIndices,nested=False):
    #Get the covering map for all patches
    if not hp.isnsideok(nside):
        raise ValueError('Incorrect map nside {0}'.format(nside))
    MapCover=np.zeros((hp.nside2npix(nside)))
    npatches=PatchIndices.shape[0]

    PatchesPosition=get2DPatchesPositions(nside,PatchWidth,PatchIndices,
                                                                nested=False)
    for kp in range(npatches):
        MapCover.flat[PatchesPosition[kp]]+=np.ones(PatchWidth*PatchWidth)
    return MapCover

def gather2DPatches(nside,Patches,PatchIndices,nested=False):
    #Get all pixel position on the sphere for all patches
    if not hp.isnsideok(nside):
        raise ValueError('Incorrect map nside {0}'.format(nside))
    PatchWidth=np.int(np.sqrt(Patches[0].size))
    MapRecons=np.zeros((hp.nside2npix(nside)))
    MapCover=np.zeros((hp.nside2npix(nside)))
    npatches=PatchIndices.shape[0]

    PatchesPosition=get2DPatchesPositions(nside,PatchWidth,PatchIndices,
                                                                nested=False)
    for kp in range(npatches):
        MapCover.flat[PatchesPosition[kp]]+=np.ones(PatchWidth*PatchWidth)
        MapRecons.flat[PatchesPosition[kp]]+=Patches[kp].flatten()
    lstCover=np.where(MapCover>0)[0]
    MapRecons[lstCover]/= MapCover[lstCover]
    return MapRecons,MapCover

def getRandomPatches(nside,ntrain,Patches,seed=None,Mask=None,nested=False):
    #Get a random subset of all patches
    np.random.seed(seed=seed)
    npatches=np.shape(Patches)[0]
    trainIndices=np.random.randint(0,high=npatches,size=ntrain,dtype='int64')
    PatchWidth=np.int(np.sqrt(Patches[0].size))

    if(Mask is not None):
        MaskPatches,ppf=getAll2DPatches(Mask,PatchWidth,nested=nested)
        npatches=MaskPatches.shape[0]
        MaskPatchesR=np.reshape(MaskPatches,(npatches,PatchWidth*PatchWidth))

        lst_rem=np.unique(np.where(MaskPatchesR[trainIndices]==0)[0])
        nrem=np.size(lst_rem)
        while(nrem>0):
            trainRemIndices=np.random.randint(0,high=npatches,size=nrem,dtype='int64')
            trainIndices[lst_rem]=trainRemIndices
            lst_rem=np.unique(np.where(MaskPatchesR[trainIndices]==0)[0])
            nrem=np.size(lst_rem)

    return Patches[trainIndices],trainIndices


def getFace2DPatches(Imag,PatchWidth,face,nested=False,BorderCap=False,
                                                                Verbose=False):
    #Get all patches from a given face
    npix=np.shape(Imag)[0]
    if not hp.isnpixok(npix):
        raise ValueError('Incorrect map size {0}'.format(npix))
    nside=hp.npix2nside(npix)
    if(nside<PatchWidth):
        raise ValueError('Incorrect PatchWidth {0} vs nside {1}'.format(
                                                       PatchWidth,nside))
    Patches2D=np.zeros((npix//12,PatchWidth,PatchWidth))
    NPatchPerSide=nside
    #Process all patches without -1
    if(Verbose):
        print("PROCESS FACE {0}".format(face))
    FaceIndices=constructExtendedFaceIndices(nside,face,PatchWidth,nested=nested)
    ExtendedFace=Imag[FaceIndices.astype('int64')]
    FacePatches2D=np.array([ExtendedFace[y:y+PatchWidth,x:x+PatchWidth]
                    for y in range(nside) for x in range(nside)
                    if -1 not in FaceIndices[y:y+PatchWidth,x:x+PatchWidth]])
    NPatchPerFace=np.shape(FacePatches2D)[0]
    return np.array(FacePatches2D),NPatchPerFace

def putFace2DPatches(Imag,Cover,Patches,face,nested=False,BorderCap=False,
                                                                Verbose=False):
    #Put all patches to a given face
    npix=np.shape(Imag)[0]
    PatchWidth=np.int64(np.sqrt(Patches.shape[1]))
    if not hp.isnpixok(npix):
        raise ValueError('Incorrect map size {0}'.format(npix))
    nside=hp.npix2nside(npix)
    if(nside<PatchWidth):
        raise ValueError('Incorrect PatchWidth {0} vs nside {1}'.format(
                                                            PatchWidth,nside))
    NPatchPerSide=nside
    #Process all patches without -1
    if(Verbose):
        print("PROCESS FACE {0}".format(face))
    FaceIndices=constructExtendedFaceIndices(nside,face,PatchWidth,nested=nested)
    FacePatches2D=np.array([np.ndarray.flatten(FaceIndices[y:y+PatchWidth,x:x+PatchWidth])
                    for y in range(nside) for x in range(nside)
                    if -1 not in FaceIndices[y:y+PatchWidth,x:x+PatchWidth]])
    for k in range(FacePatches2D.shape[0]):
        Imag[FacePatches2D[k]]=Imag[FacePatches2D[k]]+ Patches[k].flatten()
        Cover[FacePatches2D[k]]= Cover[FacePatches2D[k]]+ 1


def reconsMultiScale(LstfileName,LmaxTab,Nside=None):
    #From a list of fits filenames ordered by scales (0=highest frequency)
    #reconstruct an image, according to a maximal multipole in LmaxTab, and
    #Nside for the final map
    mapSum=hp.read_map(LstfileName[0])
    if Nside is None:
        npix=mapSum.shape[0]
        Nside=hp.npix2nside(npix)
    else:
        npix=hp.nside2npix(Nside)
    # Check that the Nside of the first map matches the desired Nside
    #If not, perform interpolation in  spherical harmonics
    if(hp.npix2nside(mapScale.shape[0]) !=Nside):
        almSum=hp.map2alm(mapScale, lmax=LmaxTab[0],mmax=LmaxTab[0])
        mapSum =hp.alm2map(almSum,Nside)

    #Set each scale to the desired NSIDE and add to the sum
    for k in range(1,len(LmaxTab)):
        mapScale=hp.read_map(LstfileName[k])
        almScale=hp.map2alm(mapScale, lmax= LmaxTab[k],mmax=LmaxTab[k])
        mapSum = mapSum + hp.alm2map(almScale,Nside)
    return mapSum

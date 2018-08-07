import numpy as np
import healpy as hp
import math

def get_all_regions():
    return get_region(lst_id=range(22))

def get_region(lst_id=None,name=None):
    #dictionary of parameters for different galactic regions used in Planck
    ref_region = ["Polaris","Orion","Pipe","Ophiuchus","Taurus","RCrA",
      "Chamaeleon-South","Pyxis", "Aquila","Auriga","RCrA-Tail","Hercules",
      "Libra","Chamaeleon-Musca","Aquila-Rift","Ara","Pisces", "Microscopium",
      "Triangulum","Perseus","Pavo","Galactic-Center"]
    lc=np.array([120.,211.,0.,354.,173.,10.,315.0,240.0,42.,145.,25.,40.,350.,
         300., 18.,336.,133.,15.,325.,143.,336.,0.])
    bc=np.array([27.,-16.,4.5,15.,-15.,-22.,-22.,12.,-15.,0.,-22.,45.,40.,-13.,
         24.,-14.,-37.,-40.,-14,-25,-28,0.])
    dl = np.array([12.,12.,5.5,12.,12.,15.,12.,25.,10.,50.,15.,15.,30.,12.,25.,
         12.,12.,12.,10.,12.,12.,12.])
    db = np.array([12.,12.,5.5,12.,12.,17.,12.,15.,10.,30.,17.,50.,30.,12.,30.,
         12.,12.,12.,7.,12.,12.,12.])

    lst_match=[]
    if(name is not None):
        lst_match=[i for i,k in enumerate(ref_region) if name in k]
        if len(lst_match) ==0:
            clm=difflib.get_close_matches(name,ref_region)
            lst_match=[i for i,k in enumerate(ref_region) if k in clm]

    if len(lst_match) ==0:
        if(lst_id is not None):
            if not hasattr(lst_id, "__len__"):
                lst_id=[lst_id]
            lst_match=[k for k in lst_id if k <len(ref_region) and k >= 0]
        else:
            print("STRING {0} NOT FOUND, USE GALACTIC CENTER INSTEAD".format(name))
            lst_match=[21]
    return [{'Name':ref_region[k],'lc':lc[k],'bc':bc[k],'dl':dl[k],
                                                'db':db[k]} for k in lst_match]

def get_pixels_in_region_new(nside, lst_id=None,region=None):
    #Get a list of pixels belonging to a given region
    try:
        dict_match= get_region(lst_id=lst_id,name=region)
        lst_coll=[]
        for kreg in range(len(dict_match)):
            print(kreg)
            lc_min=(dict_match[kreg]["lc"]-dict_match[kreg]["dl"]/2)
            lc_max=(dict_match[kreg]["lc"]+dict_match[kreg]["dl"]/2)
            bc_max=(90-(dict_match[kreg]["bc"]-dict_match[kreg]["db"]/2))*math.pi/180
            bc_min=(90-(dict_match[kreg]["bc"]+dict_match[kreg]["db"]/2))*math.pi/180
            lst_ext=hp.query_strip(nside, bc_min, bc_max, inclusive=True)
            lstang=np.asarray(hp.pix2ang(nside, lst_ext,lonlat=True))
            if(lc_min<0):
                lst_rest=np.argwhere(
                                ((lstang[0]>=0)&(lstang[0]<=lc_max))|
                                (lstang[0]>lc_min+360)&(lstang[0]<=360)
                                ).flatten()
            if(lc_max>360):
                lst_rest=np.argwhere(
                                ((lstang[0]>=360)&(lstang[0]+360<=lc_max))|
                                (lstang[0]>lc_min)&(lstang[0]<=360)
                                ).flatten()
            else:
                lst_rest=np.argwhere((lstang[0]>=lc_min)&(lstang[0]<=lc_max)).flatten()
            lstang=lstang[:,lst_rest]
            lst_coll.append(lst_ext[lst_rest])
    except Exception as inst:
        print("FAILURE: ",inst)
    return lst_coll


def extract_region(hmap,lst_id=None,region=None,keep_reso=True,minpix=128,
                              iter=0,ortho=False,fwhm_deg=None,reorder_ns=None):
    #Extract a map for the current region
    try:
        dict_match= get_region(lst_id=lst_id,name=region)
        nside=hp.get_nside(hmap)
        nmaps=hp.maptype(hmap)
        print(nmaps)
        #case single map, write it as a list
        if nmaps==0:
            pmap=[hmap]
            nmaps=1
        else:
            pmap=hmap

        if fwhm_deg is not None:
            if hasattr(fwhm_deg, "__len__"):
                assert(len(fwhm_deg) == nmaps),"Not enough fwhm: {0} vs {1}"\
                                    "map(s).".format(len(fwhm_deg),nmaps)
                fwhm_use= fwhm_deg
            else:
                fwhm_use=np.zeros([nmaps])+fwhm_deg
            lmax=np.zeros((nmaps))
            for kmap in range(nmaps):
                fl=hp.sphtfunc.gauss_beam(fwhm_use[kmap]*DtoR,lmax=2048,
                                                                     pol=False)
                llist=(np.where(fl < 1e-6))[0]
                if len(llist)==0:
                    lmax[kmap]=3*nside-1
                else:
                    lmax[kmap]=llist[0]
                if (reorder_ns is not None) and (hp.isnsideok(reorder_ns)):
                    nside=reorder_ns
                sm_map=[smooth_map_reorder_alm(pmap[kmap], fwhm_deg=fwhm_use[kmap],
                              lmax=lmax[kmap],nside_out=nside,iter=iter)\
                              for kmap in range(nmaps)]
        else:
            sm_map=pmap

        patch=[]
        for kreg in range(len(dict_match)):
            rotation=(dict_match[kreg]["lc"],dict_match[kreg]["bc"],0.)
            if keep_reso is True:
                reso_arcmin=hp.nside2resol(nside,arcmin=True)
                nxpix=np.int(max([np.ceil(dict_match[kreg]["dl"]*60./
                                                        reso_arcmin),minpix]))
                nypix=np.int(max([np.ceil(dict_match[kreg]["db"]*60./
                                                        reso_arcmin),minpix]))
                print("STAT=",reso_arcmin,nxpix,nypix,rotation)

            else:
                reso=np.min([dict_match[kreg]["dl"]/minpix,
                                                dict_match[kreg]["db"]/minpix])
                reso_arcmin=reso*60.
                nxpix=np.int(dict_match[kreg]["dl"]/reso)
                nypix=np.int(dict_match[kreg]["db"]/reso)
            a=plt.figure()
            nsub=np.int(np.max((np.ceil(np.sqrt(nmaps)),1)))
            patchreg = [np.array(hp.visufunc.gnomview(map=sm_map[kmap],
                    rot=rotation,fig=a,coord='G', xsize=nxpix,ysize=nypix,
                    reso=reso_arcmin,gal_cut=0,title=dict_match[kreg]["Name"]+
                    " map "+str(kmap),flip='astro',format='%.3g',cbar=True,
                    hold=False,margins=None,notext=False,sub=(nsub,nsub,kmap+1),
                    return_projected_map=True)) for kmap in range(nmaps)]

            patch.append(patchreg)
    except Exception as inst:
        print("FAILURE: ",inst)

    return patch

def maptoL(cmap,lmax):
    #Return a band limited version of a map
    return hp.alm2map(hp.map2alm(cmap,lmax=lmax,use_weights=True),hp.get_nside(cmap))

/******************************************************************************
 **                   Copyright (C) 2016 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Florent Sureau [rewritten from S. Beckouche Code]
 **
 **    Date:  06-07/2016
 **
 **    File:  PatchExtractor.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces classes for Patch Extraction.
 **    -----------
 **  - 10/2016: Vectorized 3D Patch Extractor (A.W.)
 **
 **
 ******************************************************************************/
 

#ifndef PatchExtractor_H
#define PatchExtractor_H
#include "TempArray.h"
#include "MatrixOper.h"
#include "IM_IO.h"
#include "gslRNG.h"

class PatchExtractor {
protected:
    size_t m_Nx;
    size_t m_Ny;
    size_t m_Nz;
    size_t m_npatches;
    size_t m_patchsize;
    dblarray* m_img;
    dblarray* m_patches;
    size_t* m_patchLoc;
    gslRNG* m_RNG;
    bool m_verbose = false;
public:

    //Constructors/Destructors
    PatchExtractor(dblarray& Img): m_patches(NULL), m_patchLoc(NULL), m_RNG(NULL) {
        setImg(Img);
    }
    virtual ~PatchExtractor() {
        if(m_patches!=NULL) delete m_patches;
        if(m_patchLoc!=NULL) free(m_patchLoc);
    }
    virtual void print(bool self = true);

    //Setters
    void setRNG(gslRNG* RNG) {m_RNG = RNG;}
    void setImg(dblarray& Img);
    void setVerbose(bool verbose) {m_verbose = verbose;}

    //Getters
    virtual dblarray* getPatches() const {return m_patches;}
    virtual std::string getSignature() const =0;
    size_t getPatchLoc(size_t PatchNb) const{
        if(PatchNb > m_npatches){
            std::cout<<"PatchExtractor:Patch Nb should be <"<<m_npatches <<
                       ", not "<<PatchNb <<std::endl;
            exit(EXIT_FAILURE);
        }
        return m_patchLoc[PatchNb];
    }

    //Processing functions
    virtual void extractPatches(size_t PatchSize,size_t Npatches)=0;

    //IOs
    void writePatches(char *FitsFileName){
        fits_write_dblarr(FitsFileName, *m_patches);}

};


class PatchExtractor2D: public PatchExtractor {
   public:
	  //Constructors/Destructors
      PatchExtractor2D(dblarray &Img):PatchExtractor(Img) {}
      virtual ~PatchExtractor2D() {}
	  //Processing functions
	  virtual void extractPatches(size_t PatchSize,size_t Npatches);
	  //Getters
	  virtual std::string getSignature() const {
                                          return std::string("PatchExtractor2D");}
	  size_t getXPatchLoc(size_t PatchNb) const{
		 if(PatchNb > m_npatches){
			std::cout<<"PatchExtractor:Patch Nb should be <"<<m_npatches <<
					 ", not "<<PatchNb <<std::endl;
			exit(EXIT_FAILURE);
		 }
         return m_patchLoc[PatchNb]%m_Nx;
	  }
	  size_t getYPatchLoc(size_t PatchNb) const{
		 if(PatchNb > m_npatches){
			std::cout<<"PatchExtractor:Patch Nb should be <"<<m_npatches <<
					 ", not "<<PatchNb <<std::endl;
			exit(EXIT_FAILURE);
		 }
         return (size_t) (m_patchLoc[PatchNb]/m_Nx);
	  }
   protected:
   private:

};

class PatchExtractor3D: public PatchExtractor {
protected:
    bool m_allChannels = false;
    bool m_multiChannelPatches = false;
    size_t m_patchSize_z;

public:
	//Constructors/Destructors
	PatchExtractor3D(dblarray &Img):PatchExtractor(Img) {}
	virtual ~PatchExtractor3D() {}
    virtual void print();

    //Processing functions
    virtual void extractPatches(size_t PatchSize,size_t Npatches);

	//Getters
    virtual std::string getSignature() const;
	size_t getXPatchLoc(size_t iPatch) const;
	size_t getYPatchLoc(size_t iPatch) const;
	size_t getZPatchLoc(size_t iPatch) const;

    //Setters
    void setAllChannels(bool allChannels);
    void setMultiChannelPatches(bool multiChannelPatches);
};
#endif

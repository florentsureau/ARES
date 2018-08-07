/*============================================================================*/
/* File Description                                                           */
/*============================================================================*/
/**
 * @file        OnlinePatchExtractor.hpp
 * @version     $Revision: 1.0 $
 * @author      A. WOISELLE
 */
/*============================================================================*/

#ifndef _OnlinePatchExtractor_HPP
#define _OnlinePatchExtractor_HPP
#include <vector>
#include <cstddef>
#include "TempArray.h"
#include "MatrixOper.h"
#include "IM_IO.h"
#include "gslRNG.h"
//#include "PatchExtractor.h"

namespace ISAP {

class OnlinePatchExtractor_Settings {
public:
    std::vector<std::string> imageList; //!< List of images
    size_t patchSize;                   //!< Spatial 2D square patch size
    bool useMultipleFiles;              //!< Sample patches from all files, or from a single random file
    size_t Nx;                          //!< Image size
    size_t Ny;                          //!< Image size
    size_t Nz;                          //!< Image size

    virtual void print() const;
    OnlinePatchExtractor_Settings(const std::vector<string>& i_imageList, size_t i_patchSize = 0,
                                  bool i_useMultipleFiles = true, size_t i_Nx = 0, size_t i_Ny = 0, size_t i_Nz = 0);
};

class OnlinePatchExtractor {
protected:
    OnlinePatchExtractor_Settings m_Settings;   //!< Settings class
protected:
    bool m_verbose = false;                     //!<
    size_t m_nFiles = 0;                        //!< Number of image files
    size_t m_nPatches = 0;                      //!< Number of patches extracted
    size_t* nPatchesPerImage;                   //!<
    int m_patchOrder;                           //!< Type of patch ordering

    dblarray* m_patches = nullptr;              //!< (nx=patchLength, ny=nPatches)
    dblarray* m_patchesMultiChannel = nullptr;  //!< (nx=nChannels, ny=2dPatchLength, nz=nPatches)
    dblarray* m_outputData = nullptr;           //!< Temporary output image data
    intarray* m_patchNorm = nullptr;            //!< Patch norm
    to_array<size_t,true>* m_patchLoc = nullptr;//!< Patch location

    gslRNG* m_RNG = nullptr;                    //!<

public:
    // Constructors/Destructors
    OnlinePatchExtractor(OnlinePatchExtractor_Settings &settings); //std::vector<std::string>& imageList, bool useMultipleFiles=false);
    virtual ~OnlinePatchExtractor();
    virtual void print() const;

    //Getters
    virtual std::string getSignature() const =0;
    virtual to_array<size_t,true>* getPatchLoc() const;
    virtual size_t getPatchLoc(size_t PatchNb) const;
    virtual inline size_t nx() const {return m_Settings.Nx;}
    virtual inline size_t ny() const {return m_Settings.Ny;}
    virtual inline size_t nz() const {return m_Settings.Nz;}
    virtual inline size_t patchSize() const {return m_Settings.patchSize;}
    gslRNG* getRNG() {return m_RNG;}

    //Setters
    void setRNG(gslRNG* RNG) {m_RNG=RNG;}
    void setVerbose(bool verbose) {m_verbose = verbose;}

    //Processing functions
    virtual dblarray* getPatchPointer(bool multiChannelPatches = false) const;
    virtual dblarray* getPatches(size_t Npatches, bool multiChannelPatches) = 0;

    //IOs
    dblarray* readFile(int iFile) const;
    dblarray* readFitsFile(int iFile) const;
    void writePatches(char *FitsFileName) const;

};


class OnlinePatchExtractor3D_Settings : public OnlinePatchExtractor_Settings {
public:
    bool allChannels = false;                 //!< if true : patchsize_z = m_Nz
    bool multiChannelPatches = false;         //!<

    virtual void print() const;
    OnlinePatchExtractor3D_Settings(const std::vector<std::string> &i_imageList, size_t i_patchSize=0,
                                    bool i_useMultipleFiles=true, bool i_allChannels = false, bool i_multiChannelPatches = false,
                                    size_t i_Nx=0, size_t i_Ny=0, size_t i_Nz=0);
};

class OnlinePatchExtractor3D: public OnlinePatchExtractor {
protected:
    OnlinePatchExtractor3D_Settings m_Settings; //!< Settings class
protected:
    size_t m_patchSize_z = 0;                   //!< Spatial/spectral : third dimension patch size
    size_t m_patchLength = 0;                   //!<
    size_t m_totalNumberOfPatches = 0;          //!< Total number of patches in one image
    size_t m_nExtractedPatches = 0;
    size_t m_currentFileNumber = -1;
    dblarray *m_centerValues = nullptr;

public:
    //Constructors/Destructors
    OnlinePatchExtractor3D(OnlinePatchExtractor3D_Settings& settings);
    ~OnlinePatchExtractor3D();
    void print() const;

    //Getters
    std::string getSignature() const;   //!< Object name string
    inline size_t nx() const {return m_Settings.Nx;}
    inline size_t ny() const {return m_Settings.Ny;}
    inline size_t nz() const {return m_Settings.Nz;}
    inline size_t patchSize() const {return m_Settings.patchSize;}
    size_t getXPatchLoc(size_t iPatch) const;   //!< patch X coordinate
    size_t getYPatchLoc(size_t iPatch) const;   //!< patch Y coordinate
    size_t getZPatchLoc(size_t iPatch) const;   //!< patch Z coordinate
    size_t getNPatchLoc(size_t PatchNb) const;  //!< patch source image number
    size_t getPatchLength() const;              //!< patch number of pixels
    bool morePatchesAvailable() const;          //!< true if not all patches have been extracted (in Natural mode)

    //Setters
    void setPatchSize(size_t i_patchSize);                  //!< Setter for m_patchSize
    void setAllChannels(bool allChannels);                  //!< Setter for m_allChannels
    void setMultiChannelPatches(bool multiChannelPatches);  //!< Setter for m_allChannels

    // get All Patches
    dblarray* getAllPatches(bool multiChannelPatches = false);
    void putAllPatches(dblarray &i_patches, dblarray &o_data);

    // get Patches
    dblarray* getRandomPatches(size_t i_nPatches, bool multiChannelPatches = false); //!< To extract random patches from imageList
    dblarray* getPatches(size_t kFile, size_t nPatches, bool multiChannelPatches = false); //!< To extract patches in natural order from a given file
    dblarray* getPatches(size_t nPatches, bool multiChannelPatches = false); //!< To extract patches in natural order from imageList
    void setPatches(dblarray *i_Patches = nullptr); //!< To update the patches of extractor object
    int putPatches(dblarray *i_Patches = nullptr);  //!< To put the patches back into the image
    int getImage(dblarray &o_data);                 //!< To get the image reconstructed from patches

    // Normalization
    void centerPatches(dblarray *Patches = nullptr);
    void inverseCenterPatches(dblarray *Patches = nullptr);

    // Deprecated get Random Patches
    dblarray* getRandomPatchesOld(size_t i_nPatches);   //!< Legacy function for random patches extraction

protected:
    void computePatchSize();

    // getPatches
    int calcRandomPatchLoc(size_t i_nPatches);                  //!<
    int calcNaturalPatchLoc(size_t kImage, size_t i_nPatches);  //!<
    int extractChosenPatches();                                 //!< Extract patches at precalculated coordinates

    // getImage
    void patchesToData_uint8(dblarray &io_data);
};

} // ISAP
#endif

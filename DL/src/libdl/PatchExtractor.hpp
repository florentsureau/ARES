/*============================================================================*/
/* File Description                                                           */
/*============================================================================*/
/**
 * @file        PatchExtractor.hpp
 * @version     $Revision: 1.0 $
 * @author      A. WOISELLE
 */
/*============================================================================*/

#ifndef _PatchExtractor_HPP
#define _PatchExtractor_HPP
#include <vector>
#include <cstddef>
#include "TempArray.h"
#include "MatrixOper.h"
#include "IM_IO.h"
#include "gslRNG.h"
//#include "PatchExtractor.h"

namespace ISAP {


class PatchExtractor_Settings {
public:
    std::vector<std::string> imageList; //!< List of images
    size_t patchSize;                   //!< Spatial 2D square patch size
    bool useMultipleFiles;              /**!< Sample patches from all files,
                                                or from a single random file*/

    virtual void print() const;
    void init();
    virtual ~PatchExtractor_Settings() {}
    PatchExtractor_Settings();
    PatchExtractor_Settings(const std::vector<string>& i_imageList,
                        size_t i_patchSize = 0, bool i_useMultipleFiles = true);
};

class PatchExtractor {
protected:
    bool m_verbose = false;                     //!<
    gslRNG* m_RNG = nullptr;                    //!<
public:
    // Constructors/Destructors
    PatchExtractor() {}
    virtual ~PatchExtractor() {}
    virtual void print() const {}

    //Getters
    virtual std::string getSignature() const {
                                   return std::string("PatchExtractor"); }

    //Setters
    virtual void setRNG(gslRNG* RNG) {m_RNG=RNG;}
    virtual void setVerbose(bool verbose) {m_verbose = verbose;}

    //Processing functions
    virtual dblarray* getPatchPointer(bool multiChannelPatches = false) const;
    virtual dblarray* getRandomPatches(size_t i_nPatches,
                                            bool multiChannelPatches = false,
                                      to_array<size_t,true>* i_patchLoc = NULL);
    virtual dblarray* getPatches(size_t kFile, size_t nPatches,
                                              bool multiChannelPatches = false);
    virtual dblarray* getPatches(size_t nPatches,
                                              bool multiChannelPatches = false);
    virtual dblarray* getAllPatches(bool multiChannelPatches = false);
    virtual void setPatches(dblarray *i_Patches = nullptr);
    virtual int putPatches(dblarray *i_Patches = nullptr);
    virtual int getImage(dblarray &o_data);
};


class PatchExtractor3D_Settings : public PatchExtractor_Settings {
public:
    bool allChannels = false;                 //!< if true : patchsize_z = m_Nz
    bool multiChannelPatches = false;         //!<

    virtual void print() const;
    void init();
    PatchExtractor3D_Settings();
    PatchExtractor3D_Settings(const std::vector<std::string> &i_imageList,
                            size_t i_patchSize=0,bool i_useMultipleFiles=true,
                bool i_allChannels = false, bool i_multiChannelPatches = false);
};

class PatchExtractor3D: public PatchExtractor {
protected:
    PatchExtractor3D_Settings m_settings; //!< Settings class
protected:
    size_t m_nFiles = 0;                        //!< Number of image files
    size_t m_nPatches = 0;                      //!< Number of patches extracted
    size_t* nPatchesPerImage;                   //!<
    int m_patchOrder;                           //!< Type of patch ordering
    bool m_usingInputImage;                     //!< Using input image (not iomageList)
    dblarray m_inputImage;                      //!< Input image if used
    size_t m_Nx = 0;                            //!< Image size
    size_t m_Ny = 0;                            //!< Image size
    size_t m_Nz = 0;                            //!< Image size
    dblarray* m_patches = nullptr;              //!< (nx=patchLength, ny=nPatches)
    dblarray* m_patchesMultiChannel = nullptr;  //!< (nx=nChannels, ny=2dPatchLength, nz=nPatches) : pointer to the same memory as m_patches //FCS: nx<->ny
    dblarray* m_outputData = nullptr;           //!< Temporary output image data
    intarray* m_patchNorm = nullptr;            //!< Patch norm
    to_array<size_t,true>* m_patchLoc = nullptr;//!< Patch location
protected:
    size_t m_patchSize_z = 0;                   //!< Spatial/spectral : third dimension patch size
    size_t m_patchLength = 0;                   //!<
    size_t m_totalNumberOfPatches = 0;          //!< Total number of patches in one image
    size_t m_nExtractedPatches = 0;
    size_t m_currentFileNumber = -1;
    dblarray *m_centerValues = nullptr;

public:
    //Constructors/Destructors
    PatchExtractor3D(PatchExtractor3D_Settings& settings);
    PatchExtractor3D(dblarray &i_inputImage, PatchExtractor3D_Settings &settings);
    virtual ~PatchExtractor3D();
    void print() const;
    void printSelf() const;

    //Getters
    dblarray* getInputImage() {return &m_inputImage;}
    std::string getSignature() const;   //!< Object name string
    inline size_t patchSize() const {return m_settings.patchSize;}
    inline size_t patchSizeZ() const {return m_patchSize_z;}
    inline size_t patchLength() const {return m_patchLength;}
    virtual inline size_t nx() const {return m_Nx;}
    virtual inline size_t ny() const {return m_Ny;}
    virtual inline size_t nz() const {return m_Nz;}
    gslRNG* getRNG() {return m_RNG;}
    to_array<size_t,true>* getPatchLoc() const;
    size_t getPatchLoc(size_t PatchNb) const;
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
    void setPatchSizeZ(size_t _zsize){m_patchSize_z=_zsize;};  //!< Set number of channels used
    void setPatchLoc(to_array<size_t,true>* patchLoc);  //!< Setter for m_allChannels

    // get All Patches
    dblarray* getAllPatches(bool multiChannelPatches = false);
    void putAllPatches(dblarray &i_patches, dblarray &o_data);

    // get Patches
    virtual dblarray* getPatchPointer(bool multiChannelPatches = false) const;
    dblarray* getRandomPatches(size_t i_nPatches,bool multiChannelPatches=false,
                 to_array<size_t,true>* i_patchLoc = NULL); //!< To extract random patches from imageList
    dblarray* getPatches(size_t kFile,size_t nPatches,
                                                bool multiChannelPatches=false); //!< To extract patches in natural order from a given file
    dblarray* getPatches(size_t nPatches, bool multiChannelPatches = false); //!< To extract patches in natural order from imageList
    void setPatches(dblarray *i_Patches = nullptr); //!< To update the patches of extractor object
    int putPatches(dblarray *i_Patches = nullptr);  //!< To put the patches back into the image
    int getImage(dblarray &o_data);                 //!< To get the image reconstructed from patches

    // Normalization
    void centerPatches(dblarray *Patches = nullptr);
    void inverseCenterPatches(dblarray *Patches = nullptr);

    // Deprecated get Random Patches
    dblarray* getRandomPatchesOld(size_t i_nPatches);   //!< Legacy function for random patches extraction

    //IOs
    dblarray* readFile(int iFile) const;
    dblarray* readFitsFile(size_t iFile) const;
    void writePatches(char *FitsFileName) const;

protected:
    void computePatchSize();

    // getPatches
    int calcRandomPatchLoc(size_t i_nPatches, to_array<size_t,true>* i_patchLoc);//!< Creates random patch positions (m_patchLoc)
    int calcNaturalPatchLoc(size_t kImage, size_t i_nPatches);  //!<
public:
    int extractChosenPatches();                                 //!< Extract patches at precalculated coordinates

    // getImage
    void patchesToData_uint8(dblarray &io_data);
};

class MultiScalePatchExtractor3D_Settings : public PatchExtractor3D_Settings {
public:
    size_t smallPatchSize = 0;                 //!< thumbnail patch size

    virtual void print() const;
    void init();
    MultiScalePatchExtractor3D_Settings();
    MultiScalePatchExtractor3D_Settings(const std::vector<std::string>
                &i_imageList, size_t i_patchSize=0, size_t i_smallPatchSize = 0,
                bool i_useMultipleFiles=true,bool i_allChannels = false,
                                            bool i_multiChannelPatches = false);
    // getPatches
    dblarray* getPatches(size_t nPatches, bool multiChannelPatches = false); //!< To extract patches in natural order from the current image
};

class MultiScalePatchExtractor3D : public PatchExtractor {
protected:
    MultiScalePatchExtractor3D_Settings m_settings; //!< Settings class
protected:
    PatchExtractor3D *m_extractor;
    PatchExtractor3D *m_thumbnailExtractor;
    dblarray m_inputImage;
    dblarray m_thumbnailImage;
    dblarray *m_MSPatches = NULL;               // (kPixel, kPatch)
    dblarray *m_MSPatchesMultiChannel = NULL;   // (kPixel, kChannel, kPatch)
    int m_patchZoom;

public:
    //Constructors/Destructors
    MultiScalePatchExtractor3D(MultiScalePatchExtractor3D_Settings &i_settings);
    MultiScalePatchExtractor3D(dblarray &i_inputImage,
                            MultiScalePatchExtractor3D_Settings &i_settings);
    ~MultiScalePatchExtractor3D();
    void print() const;
    void printSelf() const;

    dblarray* getPatchPointer(bool multiChannelPatches) const;
    dblarray createThumbnail(dblarray &i_inputImage, size_t i_patchSize,
                                                    size_t i_smallPatchSize);
    static void convX(const dblarray &M1, const dblarray &M2, dblarray &Mout);
    static void convY(const dblarray &M1, const dblarray &M2, dblarray &Mout);
    void checkThumbnailSize();

    dblarray* mergePatches(dblarray *patches, dblarray *filteredPatches,
                                            bool multiChannelPatches = false);
    dblarray* getPatches(size_t nPatches, bool multiChannelPatches = false);
    dblarray* getRandomPatches(size_t i_nPatches,bool multiChannelPatches=false,
                                    to_array<size_t,true>* i_patchLoc = NULL); //!< To extract random patches from imageList
    void setVerbose(bool verbose) {m_extractor->setVerbose(verbose);
                                    m_thumbnailExtractor->setVerbose(verbose); }
    void setRNG(gslRNG* RNG) {m_extractor->setRNG(RNG);
                                            m_thumbnailExtractor->setRNG(RNG); }
};

} // ISAP
#endif

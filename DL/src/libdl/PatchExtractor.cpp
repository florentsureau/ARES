/*============================================================================*/
/* File Description                                                           */
/*============================================================================*/
/**
 * @file        PatchExtractor.cpp
 * @version     $Revision: 1.0 $
 * @author      A. WOISELLE
 * @brief       Implementation of PatchExtractor
 */

#include <limits>

#include "PatchExtractor.hpp"
#include "ISAP_verbose.hpp"

#include "MatrixDecomposition.h"
#include "myexit.hpp"

//#define DEBUGME_(x) x
#define DEBUGME_(x)

// ======================================
// =====   PatchExtractor   =====
// =====                            =====
// ======================================

namespace ISAP {

//@comment FCS: if possible, better to use enum class than namespace+enum
namespace Extractor {
    enum PatchOrder {
        UndefinedYet,
        Random,
        Natural,
        Full
};} // Extractor


void PatchExtractor_Settings::init() {
    patchSize = 8;
    useMultipleFiles = false;
}

PatchExtractor_Settings::PatchExtractor_Settings() {
    init();
}

PatchExtractor_Settings::PatchExtractor_Settings(const std::vector<std::string>
                     &i_imageList,size_t i_patchSize, bool i_useMultipleFiles) {
    FINFO_;
    imageList = i_imageList;
    patchSize = i_patchSize;
    useMultipleFiles = i_useMultipleFiles;
    FOFNI_;
}

void PatchExtractor_Settings::print() const {
    std::cout << " imageList           = " << std::endl;
    for (std::vector<std::string>::const_iterator it = imageList.begin();
                                                    it != imageList.end(); ++it)
        std::cout << "   " << *it << std::endl;
    std::cout << " patchSize        = " << patchSize        << std::endl;
    std::cout << " useMultipleFiles = " << useMultipleFiles << std::endl;
}

dblarray* PatchExtractor::getPatchPointer(bool multiChannelPatches) const {
    std::cout << "Error : This method (" << __PRETTY_FUNCTION__ <<
                                            ") is not implemented" << std::endl;
    FINFO_;
    MyExit::exit(-1);
    return NULL;
}

dblarray* PatchExtractor::getRandomPatches(size_t i_nPatches,
                  bool multiChannelPatches, to_array<size_t,true>* i_patchLoc) {
    std::cout << "Error : This method (" << __PRETTY_FUNCTION__ <<
                                            ") is not implemented" << std::endl;
    MyExit::exit(-1);
    return NULL;
} //!< To extract random patches from imageList

dblarray* PatchExtractor::getPatches(size_t kFile, size_t nPatches,
                            bool multiChannelPatches ) {
    std::cout << "Error : This method (" << __PRETTY_FUNCTION__ <<
                                            ") is not implemented" << std::endl;
    MyExit::exit(-1);
    return NULL;
} //!< To extract patches in natural order from a given file

dblarray* PatchExtractor::getPatches(size_t nPatches, bool multiChannelPatches){
     std::cout << "Error : This method (" << __PRETTY_FUNCTION__ <<
                                            ") is not implemented" << std::endl;
     MyExit::exit(-1);
     return NULL;
} //!< To extract patches in natural order from imageList

dblarray* PatchExtractor::getAllPatches(bool multiChannelPatches){
    std::cout << "Error : This method (" << __PRETTY_FUNCTION__ <<
                                            ") is not implemented" << std::endl;
    MyExit::exit(-1);
    return NULL;
} //!< To extract all patches in natural order from imageList

void PatchExtractor::setPatches(dblarray *i_Patches ) {
    std::cout << "Error : This method (" << __PRETTY_FUNCTION__ <<
                                            ") is not implemented" << std::endl;
    MyExit::exit(-1);
} //!< To update the patches of extractor object

int PatchExtractor::putPatches(dblarray *i_Patches ) {
    std::cout << "Error : This method (" << __PRETTY_FUNCTION__ <<
                                            ") is not implemented" << std::endl;
    MyExit::exit(-1);
    return 0;
}  //!< To put the patches back into the image

int PatchExtractor::getImage(dblarray &o_data) {
    std::cout << "Error : This method (" << __PRETTY_FUNCTION__ <<
                                            ") is not implemented" << std::endl;
    MyExit::exit(-1);
    return 0;
} //!< To get the image reconstructed from patches

// ======================================
// =====   PatchExtractor3D   =====
// ======================================

void PatchExtractor3D_Settings::init() {
    allChannels = false;
    multiChannelPatches = false;
}
PatchExtractor3D_Settings::PatchExtractor3D_Settings()
    : PatchExtractor_Settings() {
    init();
}

PatchExtractor3D_Settings::PatchExtractor3D_Settings(
                    const std::vector<string>& i_imageList, size_t i_patchSize,
                    bool i_useMultipleFiles, bool i_allChannels,
                    bool i_multiChannelPatches) : PatchExtractor_Settings(
                                  i_imageList, i_patchSize,i_useMultipleFiles) {
    FINFO_;
    allChannels = i_allChannels;
    multiChannelPatches = i_multiChannelPatches;
    FOFNI_;
}

void PatchExtractor3D_Settings::print() const {
    PatchExtractor_Settings::print();
    std::cout << " allChannels        = " << allChannels        << std::endl;
    std::cout << " multiChannelPatches= " << multiChannelPatches<< std::endl;
}

void PatchExtractor3D::printSelf() const {
    //PatchExtractor::printSelf();
    std::cout << " nPatches     = " << m_nPatches   << std::endl;
    std::cout << " patchOrder   = " << m_patchOrder << std::endl;
    std::cout << " Nx           = " << m_Nx         << std::endl;
    std::cout << " Ny           = " << m_Ny         << std::endl;
    std::cout << " Nz           = " << m_Nz         << std::endl;
    std::cout << " patches*     = " << m_patches    << std::endl;
    std::cout << " patchLoc*    = " << m_patchLoc   << std::endl;
    std::cout << " RNG*         = " << m_RNG        << std::endl;

    std::cout<<" patchSize_z          = "<<m_patchSize_z        << std::endl;
    std::cout<<" patchLength          = "<<m_patchLength        << std::endl;
    std::cout<<" totalNumberOfPatches = "<<m_totalNumberOfPatches << std::endl;
    std::cout<<" nExtractedPatches    = "<<m_nExtractedPatches  << std::endl;
}

void PatchExtractor3D::print() const {
    std::cout << "PatchExtractor3D Information" << std::endl;
    m_settings.print();
    printSelf();
}

PatchExtractor3D::PatchExtractor3D(PatchExtractor3D_Settings &settings)
    : m_settings(settings) {
    FINFO_;
    m_usingInputImage = false;
    m_inputImage = dblarray();
    m_nFiles = m_settings.imageList.size();
    nPatchesPerImage = new size_t[m_nFiles];
    m_patches = new dblarray();
    m_patchesMultiChannel = new dblarray();
    m_patchNorm = new intarray();
    m_patchLoc = new  to_array<size_t,true>();
    m_outputData = new dblarray();
    m_patchOrder = Extractor::UndefinedYet;
    if(m_nFiles==0) {
        std::cout << "Error : There must be at least one file"<<
                                   " in PatchExtractor::imageList" << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    std::string imgName = m_settings.imageList[0];
    std::string ext = imgName.substr(imgName.find_last_of('.')+1);
    std::string extfits("fits");
    if(ext.compare(extfits)!=0) {
        std::cout << "ext,extfits = " << ext << "," << extfits << std::endl;
        std::cout << "Error : PatchExtractor::imageList must contain only"<<
                                                     " fits data" << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    dblarray* Img = readFile(0);
    m_Nx = Img->nx();
    m_Ny = Img->ny();
    m_Nz = max(Img->nz(),1);
    delete Img;

    computePatchSize();
    m_centerValues = new dblarray();
    // Propagate to mother class
    //PatchExtractor::m_settings = m_settings;
    FOFNI_;
}

PatchExtractor3D::PatchExtractor3D(dblarray &i_inputImage,
                    PatchExtractor3D_Settings &settings):PatchExtractor(),
                                                        m_settings(settings) {
    FINFO_;
    m_usingInputImage = true;
    m_nFiles = 1;
    nPatchesPerImage = new size_t[m_nFiles];
    m_patches = new dblarray();
    m_patchesMultiChannel = new dblarray();
    m_patchNorm = new intarray();
    m_patchLoc = new to_array<size_t,true>();
    m_outputData = new dblarray();
    m_patchOrder = Extractor::UndefinedYet;

    m_inputImage.alloc(i_inputImage.buffer(), i_inputImage.nx(),i_inputImage.ny(),
                                                            i_inputImage.nz());
    m_Nx = m_inputImage.nx();
    m_Ny = m_inputImage.ny();
    m_Nz = max(m_inputImage.nz(),1);

    computePatchSize();
    m_centerValues = new dblarray();
    // Propagate to mother class
    //PatchExtractor::m_settings = m_settings;
    FOFNI_;
}


PatchExtractor3D::~PatchExtractor3D() {
    if(nPatchesPerImage!=NULL)      delete nPatchesPerImage;
    if(m_patches!=NULL)             delete m_patches;
    if(m_patchesMultiChannel!=NULL) delete m_patchesMultiChannel;
    if(m_patchNorm!=NULL)           delete m_patchNorm;
    if(m_patchLoc!=NULL)            delete m_patchLoc;
    if(m_outputData!=NULL)          delete m_outputData;
    if(m_centerValues!=NULL)        delete m_centerValues;
}

dblarray* PatchExtractor3D::readFitsFile(size_t iFile) const {
    dblarray* pImg = new dblarray();
    if(iFile>m_settings.imageList.size()-1) {
        std::cout << __PRETTY_FUNCTION__ << " ERROR : requesting file #" << iFile << " in a list of " << m_settings.imageList.size() << " elements"  << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    fits_read_dblarr(m_settings.imageList[iFile].c_str(), *pImg);
    return pImg;
}

dblarray* PatchExtractor3D::readFile(int iFile) const {
    dblarray* pImg;
    pImg = readFitsFile(iFile);
    return pImg;
}

void PatchExtractor3D::writePatches(char *FitsFileName) const {
    fits_write_dblarr(FitsFileName, *m_patches);
}


// ======================================
// =====   PatchExtractor3D   =====
// =====          GETTERS           =====
// ======================================

std::string PatchExtractor3D::getSignature() const {
    return std::string("PatchExtractor3D");
}

//! Output pointer to patches
dblarray* PatchExtractor3D::getPatchPointer(bool multiChannelPatches) const {
    if(multiChannelPatches)
        return m_patchesMultiChannel;
    else
        return m_patches;
}

to_array<size_t,true>* PatchExtractor3D::getPatchLoc() const {
    return m_patchLoc;
}

//Getters
size_t PatchExtractor3D::getPatchLoc(size_t iPatch) const {
    if(iPatch > m_nPatches){
        std::cout << "PatchExtractor:Patch Nb should be <" << m_nPatches << ", not " << iPatch << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    return m_patchLoc->buffer()[iPatch];
}


size_t PatchExtractor3D::getXPatchLoc(size_t iPatch) const{
    if(iPatch > m_nPatches){
        std::cout << "PatchExtractor:Patch Nb should be <" << m_nPatches << ", not " << iPatch << std::endl;
        MyExit::exit(EXIT_FAILURE);
	}
    return m_patchLoc->buffer()[iPatch] % nx();
}

size_t PatchExtractor3D::getYPatchLoc(size_t iPatch) const{
    if(iPatch > m_nPatches){
        std::cout << "PatchExtractor:Patch Nb should be <" << m_nPatches << ", not " << iPatch << std::endl;
        MyExit::exit(EXIT_FAILURE);
	}
    return static_cast<size_t> ( ( m_patchLoc->buffer()[iPatch] % (nx()*ny()) ) / nx() );
}

size_t PatchExtractor3D::getZPatchLoc(size_t iPatch) const{
    if(iPatch > m_nPatches){
        std::cout << "PatchExtractor:Patch Nb should be <" << m_nPatches << ", not " << iPatch << std::endl;
        MyExit::exit(EXIT_FAILURE);
	}
    return static_cast<size_t> ( ( m_patchLoc->buffer()[iPatch] % (nx()*ny()*nz()) ) / (nx()*ny()) );
}

size_t PatchExtractor3D::getNPatchLoc(size_t iPatch) const {
    if(iPatch > m_nPatches){
        std::cout << "PatchExtractor:Patch Nb should be <" << m_nPatches << ", not " << iPatch << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    return static_cast<size_t> ( m_patchLoc->buffer()[iPatch] / (nx()*ny()*nz()) );
}

size_t PatchExtractor3D::getPatchLength() const {
    return m_patchLength;
}

//! \brief PatchExtractor3D::morePatchesAvailable
//! \return {
//!     Random mode : always false
//!     Natural mode : false if number of extracted patches < total number of patches
//! }
bool PatchExtractor3D::morePatchesAvailable() const {
    return (m_totalNumberOfPatches > m_nExtractedPatches) && (m_patchOrder==Extractor::Natural);
}

// ======================================
// =====   PatchExtractor3D   =====
// =====          SETTERS           =====
// ======================================
void PatchExtractor3D::setAllChannels(bool allChannels) {
    m_settings.allChannels = allChannels;
    computePatchSize();
}
void PatchExtractor3D::setMultiChannelPatches(bool multiChannelPatches) {
    m_settings.multiChannelPatches = multiChannelPatches;
}
void PatchExtractor3D::setPatchSize(size_t i_patchSize) {
    m_settings.patchSize = i_patchSize;
    computePatchSize();
    // Propagate to mother class
    //PatchExtractor::m_settings = m_settings;
}
void PatchExtractor3D::setPatchLoc(to_array<size_t,true>* patchLoc) {
    m_patchLoc->alloc(patchLoc->buffer(), patchLoc->nx(), patchLoc->ny(), patchLoc->nz());
}

void PatchExtractor3D::computePatchSize() {
    FINFO_;
    if(m_settings.allChannels) m_patchSize_z = nz();
    if(m_patchSize_z==0) m_patchSize_z=min(patchSize(),nz());
    m_patchSize_z=min(m_patchSize_z,nz());//FCS EDIT:do not want channel dimension
    //same as spatial dimension; m_patchSize_z can now be set
    //else m_patchSize_z = min(patchSize(), nz());
    m_patchLength = patchSize()*patchSize()*m_patchSize_z;
    m_totalNumberOfPatches = (nx()-patchSize()+1) * (ny()-patchSize()+1) * (nz()-m_patchSize_z+1);
    DEBUG_(std::cout << "patchsize_z = " << m_patchSize_z << std::endl;);
    FOFNI_;
}

// ======================================
// =====   PatchExtractor3D   =====
// =====        getPatches          =====
// ======================================

int PatchExtractor3D::calcNaturalPatchLoc(size_t kFile, size_t i_nPatches) {
    FINFO_;
    size_t *iterPatchLoc;
    // Number of patches to extract
    size_t nPatches = std::min(m_totalNumberOfPatches-m_nExtractedPatches, i_nPatches);
    if(nPatches==0) {
        ERROR_(std::cout << "All patches (" << m_nExtractedPatches << "/" <<
               m_totalNumberOfPatches << ") have already been extracted" << std::endl;);
        m_nPatches = nPatches;
        m_patchLoc->alloc(m_nPatches, 0, 0);
        return EXIT_FAILURE;
    }
    m_nPatches = nPatches;

    // Create the patch locations
    m_patchLoc->alloc(m_nPatches, 0, 0);
    iterPatchLoc = m_patchLoc->buffer();

    // first patch localization
    size_t x0 = ( m_nExtractedPatches % (nx()-patchSize()+1) );
    size_t y0 = ( m_nExtractedPatches % ((nx()-patchSize()+1)*(ny()-patchSize()+1)) ) / (nx()-patchSize()+1);
    size_t z0 = ( m_nExtractedPatches ) / ((nx()-patchSize()+1)*(ny()-patchSize()+1));
    DEBUGME_(std::cout << " patches (" << m_nExtractedPatches << ":" << m_nExtractedPatches+nPatches-1 << ")/" << m_totalNumberOfPatches-1 << std::endl;);
    DEBUGME_(std::cout << "  xyz0 =" << x0 << "," << y0 << "," << z0 << std::endl;);

    size_t k = 0; // current local patch number
    if(k<nPatches)
        for( size_t x=x0 ; x<nx()-patchSize()+1 && k<nPatches; ++x, ++k, ++iterPatchLoc ) {
            DEBUGME_( std::cout << " Selecting patch#" << k << " image patch #" << x + y0*nx() + z0*ny()*nx() + kFile*nz()*ny()*nx() <<
                      " : xyz:n " << x << "," << y0 << "," << z0 << std::endl; );
            *iterPatchLoc = x + y0*nx() + z0*ny()*nx() + kFile*nz()*ny()*nx();
        }
    if(k<nPatches)
        for( size_t y=y0+1 ; y<ny()-patchSize()+1 ; ++y )
            for( size_t x=0 ; x<nx()-patchSize()+1 && k<nPatches ; ++x, ++k, ++iterPatchLoc ) {
                *iterPatchLoc = x + y*nx() + z0*ny()*nx() + kFile*nz()*ny()*nx();
                DEBUGME_( std::cout << " Selecting patch#" << k << " image patch #" << x + y*nx() + z0*ny()*nx() + kFile*nz()*ny()*nx() <<
                          " : xyz:n " << x << "," << y << "," << z0 << std::endl; );
            }
    if(k<nPatches)
        for( size_t z=z0+1 ; z<nz()-m_patchSize_z+1 ; ++z )
            for( size_t y=0 ; y<ny()-patchSize()+1 ; ++y )
                for( size_t x=0 ; x<nx()-patchSize()+1 && k<nPatches ; ++x, ++k, ++iterPatchLoc ) {
                    *iterPatchLoc = x + y*nx() + z*ny()*nx() + kFile*nz()*ny()*nx();
                    DEBUGME_( std::cout << " Selecting patch#" << k << " image patch #" << x + y*nx() + z*ny()*nx() + kFile*nz()*ny()*nx() <<
                              " : xyz:n " << x << "," << y << "," << z << std::endl; );
                }
    DEBUGME_(m_patchLoc->display(16,16,3));
    FOFNI_;
    return EXIT_SUCCESS;
}

int PatchExtractor3D::calcRandomPatchLoc(size_t i_nPatches, to_array<size_t,true>* i_patchLoc) {
    FINFO_;
    size_t *iterPatchLoc;
    m_nPatches = i_nPatches;
    if(m_RNG==NULL) {
        ERROR_( std::cerr<<"PatchExtractor3D::calcRandomPatchLoc\n RNG not set"<<std::endl; );
        return EXIT_FAILURE;
    }
    // Number of patches
    for(size_t k=0;k<m_nFiles;k++)
        nPatchesPerImage[k] = 0;
    if(m_settings.useMultipleFiles) {
        for(size_t k=0;k<m_nPatches;k++)
            (nPatchesPerImage[m_RNG->getUniformInt(m_nFiles)])++;
    } else {
        int k = m_RNG->getUniformInt(m_nFiles);
        nPatchesPerImage[k] = m_nPatches;
    }
    // Create the random patch locations or load input positions
    if(i_patchLoc == NULL) {
        m_patchLoc->alloc(m_nPatches, 0, 0);
        iterPatchLoc = m_patchLoc->buffer();
        size_t x,y,z;
        for(size_t kFile=0 ; kFile<m_nFiles ; kFile++) {
            size_t nPatchesImg = nPatchesPerImage[kFile];
            if (nPatchesImg>0) {
                for( size_t p=0 ; p<nPatchesImg ; ++p, ++iterPatchLoc ) {
                    x = m_RNG->getUniformInt(nx()-patchSize()+1);
                    y = m_RNG->getUniformInt(ny()-patchSize()+1);
                    z = m_RNG->getUniformInt(nz()-m_patchSize_z+1);
                    *iterPatchLoc = x + y*nx() + z*ny()*nx() + kFile*nz()*ny()*nx(); // PatchLoc[kPatch] = #(x,y,z,kFile)
                } } }
    } else { // user defined positions
        if((size_t)i_patchLoc->nx() != m_nPatches) {
            std::cerr << "PatchExtractor3D::calcRandomPatchLoc" << std::endl
                      << " input patchLoc size (" << i_patchLoc->nx()
                      << ") different from nPatches(" << m_nPatches << ")" << std::endl;
            return EXIT_FAILURE;
        }
        m_patchLoc->alloc(i_patchLoc->buffer(), m_nPatches, 0, 0);
    }
    FOFNI_;
    return EXIT_SUCCESS;
}

//! Extracts random patches from all files (or a single random file) declared in the imageList or from inputImage
dblarray* PatchExtractor3D::getRandomPatches(size_t nPatches, bool multiChannelPatches, to_array<size_t,true>* i_patchLoc) {
    m_patchOrder = Extractor::Random;
    calcRandomPatchLoc(nPatches, i_patchLoc);
    extractChosenPatches();
    return getPatchPointer(multiChannelPatches);
}

//! Extracts patches in natural order from the single file declared in the imageList or from inputImage
dblarray* PatchExtractor3D::getPatches(size_t nPatches, bool multiChannelPatches) {
    if(m_nFiles !=1 ) { std::cerr << "Error : There must be only one file in fileList to use the following function" << std::endl; FINFO_;}
    return getPatches(0, nPatches, multiChannelPatches);
}

//! Extracts patches in natural order from a given file declared in the imageList or from inputImage
dblarray* PatchExtractor3D::getPatches(size_t kFile, size_t nPatches, bool multiChannelPatches) {
    // If other extraction mode or change of file : reset Natural patch parameters
    if(m_patchOrder != Extractor::Natural || m_currentFileNumber != kFile) {
        m_patchOrder = Extractor::Natural;
        m_currentFileNumber = kFile;
        m_nExtractedPatches = 0;
        m_outputData->reform(nx(), ny(), nz());
        m_outputData->init(0.f);
        m_patchNorm->reform(nx(), ny(), nz());
        m_patchNorm->init(0);
    }
    calcNaturalPatchLoc(kFile, nPatches);
    extractChosenPatches();
    m_nExtractedPatches += m_nPatches;
    return getPatchPointer(multiChannelPatches);
}

//! Patch extraction function
int PatchExtractor3D::extractChosenPatches() {
    FINFO_;
    size_t x, y, z, n;

    // Patches allocation
    m_patches->alloc(patchSize()*patchSize()*m_patchSize_z, m_nPatches, 0);
    m_patches->printInfo("m_patches");


    if(!m_settings.multiChannelPatches) // classical vectorized patches : Patches_ij = Patches(x=jPix, y=iPatch) = Patches[j][i]
        m_patchesMultiChannel->alloc(m_patches->buffer(),
                         patchSize()*patchSize()*m_patchSize_z, m_nPatches, 0);
    else // multichannel patches : Patches_ijk = Patches(x=jPix, y=iChannel, z=kPatch) = Patches[k][j][i]
        m_patchesMultiChannel->alloc(m_patches->buffer(), patchSize()*patchSize(), m_patchSize_z, m_nPatches); //! @remark so it's the same buffer as not multichannel
    DEBUGME_(m_patches->printInfo(" m_patches"););
    DEBUGME_(m_patchesMultiChannel->printInfo(" m_patchesMultiChannel"););

    // Patches extraction
    size_t nFile = std::numeric_limits<size_t>::max();
    double *iterPatch = m_patches->buffer();
    dblarray* Img = nullptr;
    double *imgBuf;
    double *iterImgBuf;

    // Reading patches, where (x,y,z) is used as the top left pixel of each patch
    for( size_t kPatch=0 ; kPatch<m_nPatches ; ++kPatch ) {
        x = getXPatchLoc(kPatch);
        y = getYPatchLoc(kPatch);
        z = getZPatchLoc(kPatch);
        n = getNPatchLoc(kPatch);
        DEBUGME_( std::cout << " Extracting " << x << "," << y << "," << z << "," << n << ",nfile=" << nFile << std::endl; );
        //m_inputImage.display(16,16,1,"m_inputImage");
        if(n != nFile) {
            //if(Img != nullptr) { std::cout << "deleting Img" << std::endl; delete Img; Img = nullptr; }
            if(m_usingInputImage) {
                DEBUGME_(std::cout << " loading inputImage : nx=" <<  m_Nx <<","<< m_Ny <<","<< m_Nz << "@" << m_inputImage.buffer() << std::endl;);
                Img = new dblarray();
                Img->alloc(m_inputImage.buffer(), m_Nx, m_Ny, m_Nz);
            } else {
                Img = readFile(n);
                DEBUGME_(std::cout << " loading image " << n << std::endl;);
            }
            imgBuf = Img->buffer();
            nFile = n;
        }
        for( size_t kz=z ; kz<z+m_patchSize_z ; ++kz )
            for( size_t ky=y ; ky<y+patchSize() ; ++ky) {
                iterImgBuf = imgBuf + x + ky*nx() + kz*ny()*nx(); // &m_img[x,ky,kz]
                for( size_t kx=x ; kx<x+patchSize() ; ++kx, ++iterImgBuf, ++iterPatch )
                    *iterPatch = *iterImgBuf; // *m_patches(kPatch) = m_img[kx,ky,kz]
            }
    }
    delete Img;
    FOFNI_;
    return EXIT_SUCCESS;
}

int PatchExtractor3D::getImage(dblarray &o_data) {
    FINFO_;
    if(m_patchOrder != Extractor::Natural) {
        std::cout << "PatchExtractor:getImage : patchOrder should be Natural(" << Extractor::Natural
                  << "), it is " << m_patchOrder << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    if(morePatchesAvailable())
        std::cout << "Warning : all patches have not been extracted" << std::endl;

    // Normalization and copy
    o_data.resize(nx(), ny(), nz());
    o_data.init(0.f);
    double *ptrData_out = o_data.buffer();
    double *ptrData = m_outputData->buffer();
    int *ptrNorm = m_patchNorm->buffer();
    for( size_t z=0 ; z<nz() ; ++z )
        for( size_t y=0 ; y<ny() ; ++y )
            for( size_t x=0 ; x<nx() ; ++x, ++ptrData, ++ptrData_out, ++ptrNorm )
                if(*ptrNorm>0)
                    *ptrData_out = *ptrData / *ptrNorm;
    FOFNI_;
    return EXIT_SUCCESS;
}

void PatchExtractor3D::setPatches(dblarray *i_Patches) {
    if(i_Patches != NULL && m_patches != i_Patches) {
        if(m_patches->nx() == i_Patches->nx()
           && m_patches->ny() == i_Patches->ny()
           && std::max(m_patches->nz(),1) ==  std::max(i_Patches->nz(),1)) {
            if(m_verbose)
                std::cout << " Updating Input Patches (different from memorized ones)" << std::endl;
            DEBUG_(std::cout << " Updating Input Patches (different from memorized ones)" << std::endl; );
            *m_patches = *i_Patches;
        } else {
            ERROR_(std::cerr << "The input patches do not match the extractor's last extraction" << std::endl;);
            ERROR_(i_Patches->printInfo("input patches"); );
            ERROR_(m_patches->printInfo("extractor patches"); );
        }
    }
}

int PatchExtractor3D::putPatches(dblarray *i_Patches) {
    FINFO_;
    if(m_patchOrder != Extractor::Natural) {
        std::cout << "PatchExtractor:putPatches : patchOrder should be Natural ("
            << Extractor::Natural<< "), it is "
            << m_patchOrder << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    // Import patches if different from internal patches
    setPatches(i_Patches);
    //FCS: Why setting back the patches, and not directly using the i_Patches?

    // Patch aggregation
    int *ptrNorm0 = m_patchNorm->buffer();
    double *ptrData0 = m_outputData->buffer();
    double *ptrPatch0 = m_patches->buffer();
    int *ptrNorm;
    double *ptrData;
    double *ptrPatch;

    size_t offset;
    for( size_t k=0 ; k<m_nPatches ; ++k )
        for( size_t dz=0 ; dz<m_patchSize_z ; ++dz )
            for( size_t dy=0 ; dy<patchSize() ; ++dy ) {
                offset = getXPatchLoc(k) + (getYPatchLoc(k)+dy)*nx() + (getZPatchLoc(k)+dz)*ny()*nx();
                DEBUGME_( std::cout << " inserting patch line xy : " << getXPatchLoc(k) << "," << (getYPatchLoc(k)+dy) << std::endl; );
                ptrData = ptrData0 + offset;
                ptrNorm = ptrNorm0 + offset;
                ptrPatch = ptrPatch0 + k*m_patchLength + dy*patchSize() + dz*patchSize()*patchSize();
                for( size_t dx=0 ; dx<patchSize() ; ++dx, ptrData++, ptrNorm++, ptrPatch++ ) {
                    *ptrData += *ptrPatch;
                    *ptrNorm += 1;
                }
            }
    //DEBUGME_( m_outputData->display(100,100,100); );
    FOFNI_;
    return EXIT_SUCCESS;
}

void PatchExtractor3D::patchesToData_uint8(dblarray &io_data) {
    FINFO_;
    to_array<uint8_t,true> norm(nx(),ny(),nz());
    io_data.resize(nx(), ny(), nz());

    // Patch aggregation
    uint8_t *ptrNorm0 = norm.buffer();
    double *ptrData0 = io_data.buffer();
    double *ptrPatch0 = m_patches->buffer();
    uint8_t *ptrNorm;
    double *ptrData;
    double *ptrPatch;

    size_t offset;
    for( size_t k=0 ; k<m_nPatches ; ++k )
        for( size_t dz=0 ; dz<m_patchSize_z ; ++dz )
            for( size_t dy=0 ; dy<patchSize() ; ++dy ) {
                offset = getXPatchLoc(k) + (getYPatchLoc(k)+dy)*nx() + (getZPatchLoc(k)+dz)*ny()*nx();
                ptrData = ptrData0 + offset;
                ptrNorm = ptrNorm0 + offset;
                ptrPatch = ptrPatch0 + k*m_patchLength + dy*patchSize() + dz*patchSize()*patchSize();
                for( size_t dx=0 ; dx<patchSize() ; ++dx, ptrData++, ptrNorm++, ptrPatch++ ) {
                    *ptrData += *ptrPatch;
                    *ptrNorm += 1;
                }
            }

    // Normalization
    ptrData = ptrData0;
    ptrNorm = ptrNorm0;
    for( size_t z=0 ; z<nz() ; ++z )
        for( size_t y=0 ; y<ny() ; ++y )
            for( size_t x=0 ; x<nx() ; ++x, ++ptrData, ++ptrNorm )
                if(*ptrNorm>0) *ptrData /= *ptrNorm;
    FOFNI_;
}


// ======================================
// =====   PatchExtractor3D   =====
// =====     getRandomPatches       =====
// ======================================
/*!
 * \brief PatchExtractor3D::extractRandomPatches
 *  Legacy function for random patch extraction
 */
dblarray* PatchExtractor3D::getRandomPatchesOld(size_t i_nPatches) {
    FINFO_;
    size_t x, y, z;
    if(m_RNG==NULL) {
        ERROR_( std::cerr<<"PatchExtractor3D::RNG not set"<<std::endl; );
    }
    m_nPatches = i_nPatches;

    // Number of patches
    m_nPatches = i_nPatches;
    for(size_t k=0;k<m_nFiles;k++)
        nPatchesPerImage[k] = 0;

    if(m_settings.useMultipleFiles) {
        for(size_t k=0;k<m_nPatches;k++)
            (nPatchesPerImage[m_RNG->getUniformInt(m_nFiles)])++;
    } else {
        int k = m_RNG->getUniformInt(m_nFiles);
        nPatchesPerImage[k] = m_nPatches;
    }

    // Patches allocation
    m_patches->alloc(patchSize()*patchSize()*m_patchSize_z, m_nPatches);
    if(!m_settings.multiChannelPatches) // classical vectorized patches : Patches_ij = Patches(x=jPix, y=iPatch) = Patches[j][i]
        m_patchesMultiChannel->alloc(m_patches->buffer(), patchSize()*patchSize()*m_patchSize_z, m_nPatches, 0);
    else // multichannel patches : Patches_ijk = Patches(x=jPix, y=iChannel, z=kPatch) = Patches[k][j][i]
        m_patchesMultiChannel->alloc(m_patches->buffer(), patchSize()*patchSize(), m_patchSize_z, m_nPatches);
    DEBUG_(m_patches->printInfo(" m_patches"););
    DEBUG_(m_patchesMultiChannel->printInfo(" m_patchesMultiChannel"););

    m_patchLoc->alloc(m_nPatches, 0, 0);

    // Patches extraction
    double *iterPatch = m_patches->buffer();
    size_t *iterPatchLoc = m_patchLoc->buffer();

    for(size_t kFile=0 ; kFile<m_nFiles ; kFile++) {
        size_t nPatchesImg = nPatchesPerImage[kFile];
        if(m_verbose) std::cout << " Extracting " << nPatchesImg << " patches from image " << m_settings.imageList[kFile] << std::endl;
        if (nPatchesImg>0) {
            dblarray* Img = readFile(kFile); //! @todo update for usingInputImage case
            double *imgBuf = Img->buffer();
            double *iterImgBuf;
            // Reading patches, where (x,y,z) is used as the top left pixel of each patch
            for( size_t p=0 ; p<nPatchesImg ; ++p, ++iterPatchLoc ) {
                x = m_RNG->getUniformInt(nx()-patchSize()+1);
                y = m_RNG->getUniformInt(ny()-patchSize()+1);
                z = m_RNG->getUniformInt(nz()-m_patchSize_z+1);
                *iterPatchLoc = x + y*nx() + z*ny()*nx() + kFile*nz()*ny()*nx(); // PatchLoc[kPatch] = #(x,y,z,kFile)
                for( size_t kz=z ; kz<z+m_patchSize_z ; ++kz )
                    for( size_t ky=y ; ky<y+patchSize() ; ++ky) {
                        iterImgBuf = imgBuf + x + ky*nx() + kz*ny()*nx(); // &m_img[x,ky,kz]
                        for( size_t kx=x ; kx<x+patchSize() ; ++kx, ++iterImgBuf, ++iterPatch )
                            *iterPatch = *iterImgBuf; // *m_patches(kPatch) = m_img[kx,ky,kz]
                    }
            }
            delete Img;
        }
        //! @todo add permutation so patches are not sorted by image
        //        gsl_permutation* Perm;
        //        Perm = m_RNG->createPermutation(nPatchesImg);
        //        new_sample_ind = gsl_permutation_get(Perm, n_loop);
    }

    // Output pointer to patches
    dblarray *patchPtr = getPatchPointer();
    FOFNI_;
    return patchPtr;
}

// ======================================
// =====   PatchExtractor3D   =====
// =====      getAllPatches         =====
// ======================================

dblarray* PatchExtractor3D::getAllPatches(bool multiChannelPatches) {
    FINFO_;
    m_patchOrder = Extractor::Full;

    // Number of patches
    size_t nPatchesImg = (nx()-patchSize()+1) * (ny()-patchSize()+1) * (nz()-m_patchSize_z+1);
    for(size_t k=0;k<m_nFiles;k++)
        nPatchesPerImage[k] = nPatchesImg;
    m_nPatches = nPatchesImg * m_nFiles;

    // Patches allocation
    DEBUG_(std::cout << " m_multiChannelPatches=" << m_settings.multiChannelPatches << std::endl;);
    DEBUG_(std::cout << " m_patchSize=" << patchSize() << ", m_patchSize_z=" << m_patchSize_z << ", m_nPatches=" << m_nPatches << std::endl;);
    m_patches->alloc(patchSize()*patchSize()*m_patchSize_z, m_nPatches);
    if(!m_settings.multiChannelPatches) // classical vectorized patches : Patches_ij = Patches(x=jPix, y=iPatch) = Patches[j][i]
        m_patchesMultiChannel->alloc(m_patches->buffer(), patchSize()*patchSize()*m_patchSize_z, m_nPatches, 0);
    else // multichannel patches : Patches_ijk = Patches(x=jPix, y=iChannel, z=kPatch) = Patches[k][j][i]
        m_patchesMultiChannel->alloc(m_patches->buffer(), patchSize()*patchSize(), m_patchSize_z, m_nPatches);
    DEBUG_(m_patches->printInfo(" m_patches"););
    DEBUG_(m_patchesMultiChannel->printInfo(" m_patchesMultiChannel"););

    // Patches extraction
    double *ptrPatch0 = m_patches->buffer();
    double *ptrPatch;

    size_t kPatch = 0;
    for(size_t kFile=0 ; kFile<m_nFiles ; kFile++) {
        if(m_verbose) std::cout << " Extracting " << nPatchesImg << " patches from image " << m_settings.imageList[kFile] << std::endl;
        dblarray* Img = readFile(kFile);
        double *ptrData0 = Img->buffer();
        double *ptrData;
        for( size_t z=0 ; z<nz()-m_patchSize_z+1 ; ++z )
            for( size_t y=0 ; y<ny()-patchSize()+1 ; ++y )
                for( size_t x=0 ; x<nx()-patchSize()+1 ; ++x, ++kPatch )
                    for( size_t dz=0 ; dz<m_patchSize_z ; ++dz )
                        for( size_t dy=0 ; dy<patchSize() ; ++dy ) {
                            ptrData = ptrData0 + x + (y+dy)*nx() + (z+dz)*ny()*nx();
                            ptrPatch = ptrPatch0 + kPatch*m_patchLength + dy*patchSize() + dz*patchSize()*patchSize();
                            for( size_t dx=0 ; dx<patchSize() ; ++dx, ptrData++, ptrPatch++ )
                                *ptrPatch = *ptrData;
                        }
        delete Img;
    }

    dblarray *patchPtr = getPatchPointer(multiChannelPatches);
    FOFNI_;
    return patchPtr;
}

void PatchExtractor3D::putAllPatches(dblarray &i_patches, dblarray &o_data) {
    FINFO_;
    if(m_patchOrder != Extractor::Full) {
        std::cout << "PatchExtractor:putAllPatches : patchOrder should be Full(" << Extractor::Full
                  << "), it is " << m_patchOrder << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    // Check size
    size_t nPatches = i_patches.ny();
    size_t nPatchesTotal = (nx()-patchSize()+1) * (ny()-patchSize()+1) * (nz()-m_patchSize_z+1);
    if( nPatches != nPatchesTotal ) {
        ERROR_( std::cerr << "Number of input patches (" << nPatches
                << ") does not correspond to data size ("
                << nx() << "x" << ny() << "x" << nz()
                << ") and patch size "
                << patchSize() << "x" << patchSize() << "x" << m_patchSize_z
                << ")" << std::endl; );
        MyExit::exit(EXIT_FAILURE);
    }

    // Patch aggregation
    o_data.alloc(nx(), ny(), nz());
    o_data.init(0.);
    double *ptrData0 = o_data.buffer();
    double *ptrPatch0 = i_patches.buffer();
    double *ptrData;
    double *ptrPatch;

    size_t kPatch = 0;
    for( size_t z=0 ; z<nz()-m_patchSize_z+1 ; ++z )
        for( size_t y=0 ; y<ny()-patchSize()+1 ; ++y )
            for( size_t x=0 ; x<nx()-patchSize()+1 ; ++x, ++kPatch )
                for( size_t dz=0 ; dz<m_patchSize_z ; ++dz )
                    for( size_t dy=0 ; dy<patchSize() ; ++dy ) {
                        ptrData = ptrData0 + x + (y+dy)*nx() + (z+dz)*ny()*nx();
                        ptrPatch = ptrPatch0 + kPatch*m_patchLength + dy*patchSize() + dz*patchSize()*patchSize();
                        for( size_t dx=0 ; dx<patchSize() ; ++dx, ptrData++, ptrPatch++ )
                            *ptrData += *ptrPatch;
                    }

    // Normalization
    size_t maxN;
    double *ptrNorm;
    dblarray xNorm(nx());
    dblarray yNorm(ny());
    dblarray zNorm(nz());

    maxN = std::min(patchSize(), nx()-patchSize()+1);
    ptrNorm = xNorm.buffer();
    for( size_t x=0 ; x<nx() ; ++x, ++ptrNorm )
        *ptrNorm = 1. / std::min(std::min(1+x, nx()-x ), maxN);

    maxN = std::min(patchSize(), ny()-patchSize()+1);
    ptrNorm = yNorm.buffer();
    for( size_t y=0 ; y<ny() ; ++y, ++ptrNorm )
        *ptrNorm = 1. / std::min(std::min(1+y, ny()-y ), maxN);

    maxN = std::min(m_patchSize_z, nz()-m_patchSize_z+1);
    ptrNorm = zNorm.buffer();
    for( size_t z=0 ; z<nz() ; ++z, ++ptrNorm )
        *ptrNorm = 1. / std::min(std::min(1+z, nz()-z ), maxN);

    ptrData = ptrData0;
    for( size_t z=0 ; z<nz() ; ++z )
        for( size_t y=0 ; y<ny() ; ++y )
            for( size_t x=0 ; x<nx() ; ++x, ++ptrData )
                *ptrData *= xNorm(x) * yNorm(y) * zNorm(z);
    FOFNI_;
}

// ======================================
// =====   PatchExtractor3D   =====
// =====      CenterPatches         =====
// ======================================

void PatchExtractor3D::centerPatches(dblarray *Patches) {
    setPatches(Patches);
    dblarray* patchesMultiChannel;
    patchesMultiChannel = getPatchPointer(m_settings.multiChannelPatches);
    SapDecomposition::subtractMeanAlongDim(*patchesMultiChannel, 0, m_centerValues);
}
void PatchExtractor3D::inverseCenterPatches(dblarray *Patches) {
    setPatches(Patches);
    dblarray* patchesMultiChannel;
    patchesMultiChannel = getPatchPointer(m_settings.multiChannelPatches);
    SapDecomposition::addAlongDim(*patchesMultiChannel, 0, *m_centerValues);
}


// =====================================
// =====================================
//    MultiScale Patch Extractor 3D
// =====================================
// =====================================

void MultiScalePatchExtractor3D_Settings::init() {
    smallPatchSize = 0;
}

MultiScalePatchExtractor3D_Settings::MultiScalePatchExtractor3D_Settings()
    : PatchExtractor3D_Settings() {
    init();
}
MultiScalePatchExtractor3D_Settings::MultiScalePatchExtractor3D_Settings(const std::vector<string>& i_imageList, size_t i_patchSize, size_t i_smallPatchSize,
                                                             bool i_useMultipleFiles, bool i_allChannels, bool i_multiChannelPatches)
    : PatchExtractor3D_Settings(i_imageList, i_patchSize,
                                      i_useMultipleFiles, i_allChannels, i_multiChannelPatches) {
    FINFO_;
    smallPatchSize = i_smallPatchSize;
    FOFNI_;
}

void MultiScalePatchExtractor3D_Settings::print() const {
    PatchExtractor3D_Settings::print();
    std::cout << " smallPatchSize     = " << smallPatchSize << std::endl;
}

void MultiScalePatchExtractor3D::printSelf() const {
    //std::cout << " m_                 = " << m_           << std::endl;
}

void MultiScalePatchExtractor3D::print() const {
    std::cout << "MultiScalePatchExtractor3D Information" << std::endl;
    m_settings.print();
    printSelf();
}
//MultiScalePatchExtractor3D::MultiScalePatchExtractor3D(MultiScalePatchExtractor3D_Settings &i_settings)
//    : PatchExtractor(), m_settings(i_settings) {
//    FINFO_;
//    m_settings.print();
//    m_extractor = new PatchExtractor3D(i_settings);
//    // Propagate to mother class
//    //MultiScalePatchExtractor3D::m_settings = m_settings;
//    FOFNI_;
//}

MultiScalePatchExtractor3D::MultiScalePatchExtractor3D(MultiScalePatchExtractor3D_Settings &i_settings)
    : m_settings(i_settings) {
    FINFO_; // untested
    // Init patch extractors
    m_settings.print();
    checkThumbnailSize();
    m_extractor = new PatchExtractor3D(m_settings);
    dblarray *inputImage = m_extractor->readFile(0);
    inputImage->printInfo();
    m_thumbnailImage = createThumbnail(*inputImage, m_settings.patchSize, m_settings.smallPatchSize);
    m_thumbnailImage.printInfo("thumbnailImage");
    DEBUGME_(fits_write_dblarr("msExtractor_thumbnail.fits", m_thumbnailImage););
    m_thumbnailExtractor = new PatchExtractor3D(m_thumbnailImage, m_settings);
    m_MSPatches = new dblarray();
    m_MSPatchesMultiChannel = new dblarray();
    FOFNI_;
}

MultiScalePatchExtractor3D::MultiScalePatchExtractor3D(dblarray &i_inputImage, MultiScalePatchExtractor3D_Settings &i_settings)
    : PatchExtractor(), m_settings(i_settings) {
    FINFO_;
    m_settings.print();
    checkThumbnailSize();
    //PatchExtractor3D_Settings settings(i_settings);
    //settings.print();
    m_extractor = new PatchExtractor3D(i_inputImage, m_settings);
    m_inputImage.alloc(i_inputImage.buffer(), i_inputImage.nx(), i_inputImage.ny(), i_inputImage.nz());
    m_thumbnailImage = createThumbnail(m_inputImage, m_settings.patchSize, m_settings.smallPatchSize);
    DEBUGME_(fits_write_dblarr("msExtractor_thumbnail.fits", m_thumbnailImage););
    m_thumbnailExtractor = new PatchExtractor3D(m_thumbnailImage, m_settings);
    m_MSPatches = new dblarray();
    m_MSPatchesMultiChannel = new dblarray();
    FOFNI_;
}

MultiScalePatchExtractor3D::~MultiScalePatchExtractor3D() {
    delete m_extractor;
    delete m_thumbnailExtractor;
    delete m_MSPatches;
}

void MultiScalePatchExtractor3D::checkThumbnailSize() {
    if( m_settings.patchSize % 2 == 0 || m_settings.smallPatchSize % 2 == 0 ) {
        std::cout << "Error in " << __PRETTY_FUNCTION__ << " : " << std::endl;
        std::cout << " patchSize (" <<m_settings.patchSize<< ") and smallPatchSize (" <<m_settings.smallPatchSize<< ") must be odd" << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    double f = double(m_settings.patchSize - 1) / double(m_settings.patchSize - 1);
    if( f-floor(f) > 0 ) {
        std::cout << "Error in " << __PRETTY_FUNCTION__ << " : " << std::endl;
        std::cout << " patchSize-1 (" <<m_settings.patchSize-1<< ") and smallPatchSize-1 (" <<m_settings.smallPatchSize-1<< ") must be multiples" << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
}

dblarray MultiScalePatchExtractor3D::createThumbnail(dblarray &i_inputImage, size_t i_patchSize, size_t i_smallPatchSize) {
    dblarray thumbnail;
    dblarray temp;
    m_patchZoom = (i_patchSize-1) / (i_smallPatchSize-1) ;
    int f2 = m_patchZoom/2;
    dblarray filter = dblarray(2*f2+1);
    for(int k=0;k<2*f2+1;++k)
        filter(k) = f2+1 - abs(f2-k);
    DEBUGME_(filter.display("filter"););
    double s = 0;
    for(int k = 0 ; k < filter.n_elem(); ++k)
        s += filter(k);
    for(int k = 0 ; k < filter.n_elem(); ++k)
        filter(k) /= s;
    DEBUGME_(filter.display("filter normalized"););

    ISAP::MultiScalePatchExtractor3D::convX(i_inputImage, filter, temp);
    ISAP::MultiScalePatchExtractor3D::convY(temp, filter, thumbnail);
    return thumbnail;
}

void MultiScalePatchExtractor3D::convX(const dblarray &M1, const dblarray &M2, dblarray &Mout) {
    Mout.alloc(M1.nx(), M1.ny(), M1.nz());
    Mout.init(0.f);
    int N = M1.nx();
    int Nx2 = (M2.n_elem()-1)/2;
    int isat;
    for(int k = 0; k<max(M1.nz(),1) ; ++k)
        for(int j = 0; j<M1.ny() ; ++j)
            for(int i = 0; i<M1.nx() ; ++i)
                for(int l = -Nx2; l<=Nx2 ; ++l) {
                    isat = min(max(i-l,0),N-1);
                    Mout(i,j,k) += M1(isat,j,k) * M2(l+Nx2);
                }
}
void MultiScalePatchExtractor3D::convY(const dblarray &M1, const dblarray &M2, dblarray &Mout) {
    Mout.alloc(M1.nx(), M1.ny(), M1.nz());
    Mout.init(0.f);
    int N = M1.ny();
    int Nx2 = (M2.n_elem()-1)/2;
    int jsat;
    for(int k = 0; k<max(M1.nz(),1) ; ++k)
        for(int j = 0; j<M1.ny() ; ++j)
            for(int i = 0; i<M1.nx() ; ++i)
                for(int l = -Nx2; l<=Nx2 ; ++l) {
                    jsat = min(max(j-l,0),N-1);
                    Mout(i,j,k) += M1(i,jsat,k) * M2(l+Nx2);
                }
}

dblarray* MultiScalePatchExtractor3D::getPatchPointer(bool multiChannelPatches) const {
    if(multiChannelPatches)
        return m_MSPatchesMultiChannel;
    else
        return m_MSPatches;
}

/*!
 * \brief MultiScalePatchExtractor3D::getRandomPatches
 *     Extracts random patches at two scales and fuses them in a single vector
 * \param nPatches             : Number of patches to extract
 * \param multiChannelPatches  : To set multiChannel ordering on the output MS patches
 * \return
 */
dblarray* MultiScalePatchExtractor3D::getRandomPatches(size_t nPatches, bool multiChannelPatches, to_array<size_t,true>* i_patchLoc) {

    // extract patches
    dblarray* patches = m_extractor->getRandomPatches(nPatches, false, i_patchLoc);
    dblarray* filteredPatches = m_thumbnailExtractor->getRandomPatches(nPatches, false, m_extractor->getPatchLoc());

    //patches->display("getRandom : patches");
    //m_extractor->getPatchLoc()->display("m_extractor->getPatchLoc()");
    //m_thumbnailExtractor->getPatchLoc()->display("m_thumbnailExtractor->getPatchLoc()");

    // decimate filteredPatches and merge with patches
    //*filteredPatches *= 4.0f;
    return mergePatches(patches, filteredPatches, multiChannelPatches);
}

/*!
 * \brief MultiScalePatchExtractor3D::getPatches
 *     Extracts patches in natural order at two scales and fuses them in a single vector
 * \param nPatches             : Number of patches to extract
 * \param multiChannelPatches  : To set multiChannel ordering on the output MS patches
 * \return
 */
dblarray* MultiScalePatchExtractor3D::getPatches(size_t nPatches, bool multiChannelPatches) {
    // extract patches
    dblarray* patchesMC = m_extractor->getPatches(nPatches, false);
    dblarray* filteredPatchesMC = m_thumbnailExtractor->getPatches(nPatches, false);

    // decimate filteredPatches and merge with patches
    mergePatches(patchesMC, filteredPatchesMC, multiChannelPatches);
    return getPatchPointer(multiChannelPatches);
}

/*!
 * \brief MultiScalePatchExtractor3D::mergePatches
 *     Merges extracted patches
 * \param patches              : input high resolution patches, in non multiChannel ordering
 * \param filteredPatches      : input lower resolution (blurry, same size) patches, in non multiChannel ordering
 * \param multiChannelPatches  : To set multiChannel ordering on the output MS patches
 * \return Merged patches
 */
dblarray* MultiScalePatchExtractor3D::mergePatches(dblarray *patches, dblarray *filteredPatches, bool multiChannelPatches) {
    DEBUGME_(patches->display("Patches"););
    m_thumbnailExtractor->print();
    DEBUGME_(filteredPatches->display("thumbnailPatches"););

    // Concatenated patches
    // 2d small and large patch length
    int small2DPatchLength =  m_settings.smallPatchSize *  m_settings.smallPatchSize;
    int large2DPatchLength =  m_settings.patchSize *  m_settings.patchSize;
    std::cout << "patches sizes : " << std::endl
              << "(" << patches->nx() << "," << patches->ny() << "," << patches->nz() << ") "
              << std::endl;
    size_t nChannels = std::max(m_extractor->patchSizeZ(), static_cast<size_t>(1));
    size_t nPatches = patches->ny();
    m_MSPatchesMultiChannel->alloc(
                large2DPatchLength + small2DPatchLength,
                nChannels,
                nPatches);
    m_MSPatches->alloc(m_MSPatchesMultiChannel->buffer(),
                       ( large2DPatchLength + small2DPatchLength ) * nChannels,
                       nPatches
                       );

    // put patches in multichannel mode
    patches->reform(large2DPatchLength, nChannels, nPatches);
    filteredPatches->reform(large2DPatchLength, nChannels, nPatches);
    std::cout << "MSpatch sizes : " << std::endl
              << "(" << m_MSPatchesMultiChannel->nx() << "," << m_MSPatchesMultiChannel->ny() << "," << m_MSPatchesMultiChannel->nz() << ") "
              << "(" << m_MSPatches->nx() << "," << m_MSPatches->ny() << "," << m_MSPatches->nz() << ") "
              << std::endl;

    // Add full resolution patch
    for(size_t k = 0;  k < nPatches ; ++k)
        for(size_t j = 0; j < nChannels ; ++j)
            for(size_t i = 0; i < (size_t) large2DPatchLength ; ++i)
                m_MSPatchesMultiChannel->operator()(i,j,k) = patches->operator()(i,j,k);
    DEBUGME_(m_MSPatchesMultiChannel->display("MS Patches"););

    // Add downsampled filtered patch
    int jump = m_patchZoom;
    int i, is, i1s, i2s;
    //for(int k = 0;  k < max(patches->nz(),1) ; ++k) {
    for(size_t k = 0;  k < nPatches ; ++k) {
        for(size_t j = 0; j < nChannels ; ++j) {
            // multi-spectral offset
            i2s = 0;
            // loop on
            for(size_t i2 = 0;  i2 < m_settings.patchSize ; i2+=jump, ++i2s) {
                i1s = 0;
                for(size_t i1 = 0;  i1 < m_settings.patchSize ; i1+=jump, ++i1s) {
                    i = i1 + i2*m_settings.patchSize;
                    is = i1s + i2s*m_settings.smallPatchSize;
                    DEBUGME_(std::cout << "MSpatches(" << large2DPatchLength << "+" << is << "," << j << "," << k<< ") <- " <<
                             "filteredPatches( " << "(" << i1 << "," << i2 << ")=" << i << "," << j << "," << k << ")" << std::endl;);
                    m_MSPatchesMultiChannel->operator()(large2DPatchLength+is,j,k) = filteredPatches->operator()(i,j,k);
                }
            }
        }
    }
    DEBUGME_(m_MSPatchesMultiChannel->display("MS Patches multi channel"););

    return getPatchPointer(multiChannelPatches);
}



//    int smallContiguousPatchLength =  small2DPatchLength;
//    if(!multiChannelPatches)
//        smallContiguousPatchLength *= std::max(m_extractor->patchSizeZ(), static_cast<size_t>(1));
//    m_MSPatches = new dblarray(patches->nx() + smallContiguousPatchLength, patches->ny(), patches->nz());
// Add full resolution patch
//for(int k = 0;  k < max(patches->nz(),1) ; ++k)
//    for(int j = 0; j < nPatches ; ++j)
//        for(int i = 0; i < patches->nx() ; ++i)
//            m_MSPatchesMultiChannel->operator()(i,j,k) = patches->operator()(i,j,k);
//DEBUGME_(m_MSPatchesMultiChannel->display("MS Patches"););
//
// Add downsampled filtered patch
//    int jump = m_patchZoom;
//    int i, is, i1s, i2s;
//    if(!multiChannelPatches) {
//        for(int j = 0; j < patches->ny() ; ++j) {
//            for(size_t i3 = 0; i3 < m_extractor->patchSizeZ() ; i3++) {
//                // multi-spectral offset
//                int i0s = i3*small2DPatchLength;
//                int i0 = i3*large2DPatchLength;
//                i2s = 0;
//                // loop on
//                for(size_t i2 = 0;  i2 < m_settings.patchSize ; i2+=jump, ++i2s) {
//                    i1s = 0;
//                    for(size_t i1 = 0;  i1 < m_settings.patchSize ; i1+=jump, ++i1s) {
//                        i = i1 + i2*m_settings.patchSize;
//                        is = i1s + i2s*m_settings.smallPatchSize;
//                        DEBUGME_(std::cout << "MSpatches(" << patches->nx() << "+" << i0s << "+" << is << "," << j << "," << 0 << ") <- " <<
//                                 "filteredPatches( " << i0 << "+(" << i1 << "," << i2 << ")=" << i << "," << j << "," << 0 << ")" << std::endl;);
//                        m_MSPatchesMultiChannel->operator()(patches->nx()+i0s+is,j,0) = filteredPatches->operator()(i0+i,j,0);
//                    }
//                }
//            }
//        }
//    } else {// multispectral
//        for(int k = 0;  k < max(patches->nz(),1) ; ++k) {
//            //int i0 = k*large2DPatchLength;
//            for(int j = 0; j < patches->ny() ; ++j) {
//                    // multi-spectral offset
//                    i2s = 0;
//                    // loop on
//                    for(size_t i2 = 0;  i2 < m_settings.patchSize ; i2+=jump, ++i2s) {
//                        i1s = 0;
//                        for(size_t i1 = 0;  i1 < m_settings.patchSize ; i1+=jump, ++i1s) {
//                            i = i1 + i2*m_settings.patchSize;
//                            is = i1s + i2s*m_settings.smallPatchSize;
//                            DEBUGME_(std::cout << "MSpatches(" << patches->nx() << "+" << is << "," << j << "," << k<< ") <- " <<
//                                     "filteredPatches( " << "(" << i1 << "," << i2 << ")=" << i << "," << j << "," << k << ")" << std::endl;);
//                            m_MSPatchesMultiChannel->operator()(patches->nx()+is,j,k) = filteredPatches->operator()(i,j,k);
//                        }
//                    }
//            }
//        }
//    }



} // ISAP

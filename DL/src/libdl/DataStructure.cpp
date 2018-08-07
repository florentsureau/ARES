/******************************************************************************
 **                   Copyright (C) 2016 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Florent Sureau
 **
 **    Date:  06-07/2016
 **
 **    File:  DataStructure.cpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Implementation of data structure.
 **    -----------
 **  - 06-07/2017: added code wrapping TempArray and/or PatchExtractor3D
 ******************************************************************************/

#include "DataStructure.h"

//Overloading of operator << to print easily Manifolds
std::ostream& operator <<(std::ostream& stdout, Manifold _m){
    switch(_m){
        case Manifold::Rn   : stdout << "R";    break;
        case Manifold::Rpn: stdout << "R+"; break;
        case Manifold::S1n : stdout << "S1";  break;
        case Manifold::S1Rn  : stdout << "S1R";   break;
        case Manifold::S1wRn  : stdout << "S1wR";   break;
        default    : stdout.setstate(std::ios_base::failbit);
    }
    return stdout;
}

// ===================
//     CONSTRUCTORS
// ===================

DataStruc::DataStruc(ISAP::PatchExtractor3D* _p, ExtractorType _m, size_t npatch){
    m_extractor=_p;
    m_extract_mode=_m;
    m_nPatches=npatch;
    m_inputPatches=nullptr;
    m_seed=0;
    m_verb=false;
    m_extracted=false;
    m_contiguousChannels=false;
    m_RNG=nullptr;
    if(m_extractor!=nullptr){
        m_PatchSize=m_extractor->getPatchLength();
        m_2DPatchSize=m_extractor->patchSize();
        m_ZPatchSize=m_extractor->patchSizeZ();
        m_ownPtr=false;
    } else {
        std::cout<<"BEWARE: no extractor set : cannot extract."<<std::endl;
        m_PatchSize=0;
        m_2DPatchSize=0;
        m_ZPatchSize=0;
        m_ownPtr=false;
    }
}

DataStruc::DataStruc(dblarray* inputData,std::vector<size_t>& ChannelDim,
        std::vector<Manifold>&  ManifoldType, bool isContiguous, bool OwnPtr){
    m_extractor=nullptr;
    m_PatchSize=std::accumulate(ChannelDim.begin(), ChannelDim.end(), 0);
    m_ZPatchSize=ChannelDim.size();
    m_2DPatchSize=m_PatchSize/m_ZPatchSize;
    m_nPatches=(size_t)round(inputData->n_elem()/m_PatchSize);
    if(m_nPatches*m_PatchSize != (size_t) inputData->n_elem()){
        std::cout<<"Mismatch between input dimensions and array"<<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    m_inputPatches=new dblarray();
    m_inputPatches->alloc(inputData->buffer(),inputData->n_elem(),0,0);
    m_ownPtr=OwnPtr;
    m_verb=false;
    m_extracted=true;
    m_RNG=nullptr;
    m_seed=0;
    setManifold(ManifoldType);
    setManifoldDim(ChannelDim);
    m_contiguousChannels=isContiguous;
    if(getNManifold()==1) m_contiguousChannels=true;
}

// ===================
//     GETTERS
// ===================

string DataStruc::getSignature() const {
    std::string ClassSign("Data");
    if(m_contiguousChannels) ClassSign+=std::string("MultiPatchChannels");
    else ClassSign+=std::string("MultiChannelPatches");
    stringstream strstr;
    strstr<<ClassSign;
    for(size_t kdim=0;kdim<m_manifold.size();++kdim){
        strstr<<m_manifold[kdim];
/*        switch (m_manifold[kdim]) {
            case Manifold::Rn:
                strstr<<"R";
                break;
            case Manifold::Rpn:
                strstr<<"R+";
                break;
            case Manifold::S1n:
                strstr<<"S";
                break;
            case Manifold::S1wRn:
                strstr<<"SwR";
                break;
            default:
                exit(-1);
                break;
        }*/
        strstr<<m_productDim[kdim];
        if(kdim<m_manifold.size()-1) strstr<<"x";
    }
    if(m_extractor!=nullptr) strstr<<" with "<<m_extractor->getSignature();
    return strstr.str();
}

dblarray* DataStruc::getManifoldData(size_t Channel){
    if(Channel>=getNManifold()){
        std::cout<<"Incorrect Channel"<<Channel<<" vs "<<m_manifold.size()<<
                                                                         std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    if(!m_extracted) extractPatches();
    dblarray* ManifoldData;
    if(m_contiguousChannels){
        ManifoldData=new dblarray();
        ManifoldData->alloc(getSampleChannel(0,Channel),
                                            m_productDim[Channel],m_nPatches,0);
    } else {
        size_t Nx=getManifoldDimSize(Channel);
        size_t Ny=m_nPatches;
        ManifoldData=new dblarray(Nx,Ny);
        double* ManBuffer=ManifoldData->buffer();
        double* InternBufferStart=m_inputPatches->buffer()+m_offsetDim[Channel];
        double* InternBuffer=InternBufferStart;
        for(size_t kn=0;kn<Ny;++kn,InternBufferStart+=m_PatchSize){
            InternBuffer=InternBufferStart;
            for(size_t kp=0;kp<Nx;++kp,++InternBuffer,++ManBuffer)
                                                        *ManBuffer=*InternBuffer;
        }
    }
    return ManifoldData;
}

// ===================
//     SETTERS
// ===================


void DataStruc::resetManifold(std::vector<Manifold> _manifold, std::vector<size_t> _dim) {
    size_t msize=_manifold.size();
    if(msize!=_dim.size()){
        std::cout<<"Incorrect size inputs"<<_manifold.size()<<" vs "<<_dim.size()<<
                                                                         std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    size_t _patchSize=floor(sqrt(_dim[0]));
    if(_patchSize*_patchSize!=_dim[0]){
        std::cout<<"First dimension not square:"<<_dim[0]<< std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    if(!m_offsetDim.empty()) m_offsetDim.clear();
    m_offsetDim.push_back(0);
    for(size_t kprod=1;kprod<msize;++kprod){
        m_offsetDim.push_back(m_offsetDim[kprod-1]+_dim[kprod-1]);
    }
    //Check total dimension with Patch extractor
    size_t totalDim=m_offsetDim[msize-1]+_dim[msize-1];
    if(totalDim!=getPatchSize()){
        std::cout<<"Manifold/Patch Dimensions do not match:"<<totalDim<<" vs "<<
                                                        getPatchSize()<<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    //Check spatial dimension with Patch extractor
    if(_patchSize!=get2DPatchSize()){
        //Simple warning: we authorize the first manifold to be a cartesian product
        //of size greater than the patch size (typically n times)
        std::cout<<"Beware: spatial dimension do not match:"<<_patchSize<< "vs"<<
                                                    get2DPatchSize()<<std::endl;
    }
    if(!m_productDim.empty()) m_productDim.clear();
    if(!m_manifold.empty()) m_manifold.clear();
    m_productDim=_dim;
    m_manifold=_manifold;
}

void DataStruc::setContiguousChannels(bool ContiguousFlag){
    if((ContiguousFlag)&&(!m_contiguousChannels)){
        if(m_extracted){//Case we need to change to contiguous
            size_t Nelems=(size_t)m_inputPatches->n_elem();
            dblarray* newArray= new dblarray(Nelems);
            double* newBuffer=newArray->buffer();
            double* ItNewBufferStart=newBuffer;
            double* BufferStart=m_inputPatches->buffer();
            for(size_t kchan=0;kchan<getpatchSizeZ();++kchan){
                BufferStart=m_inputPatches->buffer()+m_offsetDim[kchan];
                for(size_t kpatch=0;kpatch<m_nPatches;kpatch++,BufferStart+=m_PatchSize){
                    double* ItBuffer=BufferStart;
                    for(size_t kpix=0;kpix<m_productDim[kchan];++kpix,++ItBuffer
                                                            ,++ItNewBufferStart)
                        *ItNewBufferStart=*ItBuffer;
                }
            }
            if(m_ownPtr) delete m_inputPatches;
            m_inputPatches=newArray;
            m_ownPtr=true;
            std::cout<<"NOW OWN PTR:"<<std::endl;
        }
    } else if((!ContiguousFlag)&&(m_contiguousChannels)){
        if(m_extracted){//Case we need to change from contiguous
            size_t Nelems=(size_t)m_inputPatches->n_elem();
            dblarray* newArray= new dblarray(Nelems);
            double* newBuffer=newArray->buffer();
            double* ItNewBufferStart=newBuffer;
            double* BufferStart=m_inputPatches->buffer();
            for(size_t kchan=0;kchan<getpatchSizeZ();++kchan){
                BufferStart=m_inputPatches->buffer()+m_offsetDim[kchan];
                for(size_t kpatch=0;kpatch<m_nPatches;kpatch++,BufferStart+=m_PatchSize){
                    double* ItBuffer=BufferStart;
                    for(size_t kpix=0;kpix<m_productDim[kchan];++kpix,++ItBuffer
                                                            ,++ItNewBufferStart)
                        *ItBuffer=*ItNewBufferStart;
                }
            }
            if(m_ownPtr) delete m_inputPatches;
            m_inputPatches=newArray;
            m_ownPtr=true;
            std::cout<<"NOW OWN PTR"<<std::endl;
        }
    }
    m_contiguousChannels=ContiguousFlag;
}

void DataStruc::putManifoldData(dblarray* Data, size_t Channel){
    if(Channel>=getNManifold()){
        std::cout<<"Incorrect Channel"<<Channel<<" vs "<<m_manifold.size()<<
                                                                         std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    if(!m_extracted){
        std::cout<<"Data not yet extracted:"<<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    size_t Nx=getManifoldDimSize(Channel);
    size_t Ny=m_nPatches;
    if(((size_t)Data->nx()!=Nx)||((size_t)Data->ny()!=Ny)){
        std::cout<<"Data to set do not match manifold size:"<<Data->nx()<<"x"<<
                                    Data->ny()<<" vs "<<Nx<<"x"<<Ny<<std::endl;
                                    MyExit::exit(EXIT_FAILURE);
    }
    double* ManBuffer=Data->buffer();
    if(m_contiguousChannels){
        double* ItBuffer=m_inputPatches->buffer()+m_offsetDim[Channel]*m_nPatches;
        for(size_t kel=0;kel<(size_t)Data->n_elem();++kel,++ManBuffer,++ItBuffer)
                                                                  *ItBuffer=*ManBuffer;
    } else {
        double* BufferStart=m_inputPatches->buffer()+m_offsetDim[Channel];
        double* ItBuffer=BufferStart;
        for(size_t kn=0;kn<Ny;++kn,BufferStart+=m_PatchSize){
            ItBuffer=BufferStart;
            for(size_t kp=0;kp<Nx;++kp,++ItBuffer,++ManBuffer) *ItBuffer=*ManBuffer;
        }
    }
}

// ===================
//     PROCESSING
// ===================

dblarray* DataStruc::extractPatches(){
    if(m_extractor==nullptr){
        std::cout<<"No extractor set : cannot extract."<<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    dblarray* Patches;
    if(m_inputPatches!=nullptr) delete m_inputPatches;

    switch (m_extract_mode) {
        case ExtractorType::Random :
            Patches = m_extractor->getRandomPatches(m_nPatches,false);
            break;
        case ExtractorType::Natural :
            Patches = m_extractor->getPatches(m_nPatches,false);
            break;
        case ExtractorType::Full :
            Patches = m_extractor->getAllPatches(false);
            break;
        default:
            exit(-1);
            break;
    }
    m_nPatches=Patches->ny();
    //Copy an array either owned or not
    if(m_contiguousChannels){
        std::cout<<"Put to Contiguous Channels"<<std::endl;
        //Extractor gives non-contiguous channels
        m_inputPatches=Patches;
        m_contiguousChannels=false;
        //Create a new array, owned with correct ordering
        setContiguousChannels(true);
    } else {
        std::cout<<"Create fake array"<<std::endl;
        m_inputPatches=new dblarray();
        //Create a new array, owned with correct structure, without modifying
        //structure of patch extractor internal array
        m_inputPatches->alloc(Patches->buffer(),Patches->n_elem(),0,0);
    }
    m_extracted=true;
    std::cout<<"DONE"<<std::endl;
    return m_inputPatches;
}

int DataStruc::gatherPatches(){
    if(m_extractor==nullptr){
        std::cout<<"No extractor set : cannot extract."<<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    //Patch Extractor has not contiguous channels
    if(m_contiguousChannels) setContiguousChannels(false);

    dblarray* Patches;
    Patches->alloc(m_inputPatches->buffer(),m_2DPatchSize*m_ZPatchSize,m_nPatches);

    if((m_extract_mode==ExtractorType::Natural)||
                                        (m_extract_mode==ExtractorType::Full))
        return m_extractor->putPatches(Patches);
    else {
        std::cout<<"Can only Gather Patches in natural or full order"<<std::endl;
        return EXIT_FAILURE;
    }
}

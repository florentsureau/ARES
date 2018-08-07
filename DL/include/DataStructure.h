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
 **    Date:  06-07/2017
 **
 **    File:  DataStructure.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces  classes for Data Structure
 **    -----------
 **  - 06-07/2017: added class to wrap TempArray and/or PatchExtractor3D
 **
 ******************************************************************************/

#ifndef DataStructure_H
#define DataStructure_H

#include "TempArray.h"
#include "MatrixDecomposition.h"
#include "IM_IO.h"
#include "ISAP_verbose.hpp"
#include "myexit.hpp"
#include "PatchExtractor.hpp"

/**
 * @enum Manifold
 * @brief choice among manifolds for data
*/
enum class Manifold {Rn,Rpn,S1n,S1Rn,S1wRn};
//R^n,R^{+n}, S1^n, S1 extended to Rn, S1 aswrapped Rn,

std::ostream& operator <<(std::ostream& stdout, Manifold _m);


enum class ExtractorType {UndefinedYet,Random, Natural,Full};

/**
 * @class DataStruc
 * @brief Class describing the structure of data considered as the product of
   different manifolds. Assumption is that a patch is either multichannel contiguous in a
   patch, i.e. for a RGB sample, we would have in one patch: all channel R, all channel G,
   all channel B and then these patches are concatenated across patches ; or the data is
   fully multichannel contiguous: we have all patches for a channel contiguous
   @see DataStruc::m_contiguousChannels flag.
*/
class DataStruc {
    protected:
        std::vector<Manifold> m_manifold;//!< Type of manifold for each channel
        std::vector<size_t> m_productDim;//!< number of points in each manifold
        std::vector<size_t> m_offsetDim;//!< number of points before each channel
        dblarray* m_inputPatches;//!< Pointer to extracted patches [nPix2D*Nchan*Npatches]
        gslRNG* m_RNG;//!< Random number generator
        size_t m_seed;//!< Seed for RNG
        size_t m_nPatches;//!< Number of patches
        size_t m_PatchSize;//!<Size of the Patch (Spatial*Channels),total dimensionality
        size_t m_2DPatchSize;//!<Spatial size of the Patch
        size_t m_ZPatchSize;//!<Number of channels in the Patch
        bool m_verb;//!< Verbosity flag
        bool m_extracted;//!< flag set if data has been extracted
        bool m_ownPtr;
        bool m_contiguousChannels=false;/**<flag to specify how data is stored;
        by default for a patch data on the same manifold is contiguous, but disjoint across
        patches (e.g. [Npix2d,NManifolds,Npatches]). If true, data on
        same manifold are contiguous (e.g. [Npix2d,Npatches,NManifolds])).
        Note however that Manifolds need not be of same size Npix2d as suggested by the
        examples.*/
    public:
        ISAP::PatchExtractor3D* m_extractor;//!< Extractor to get patches
        ExtractorType m_extract_mode;//!< Type of patch extraction
        /**
          * Constructors.
        */
        DataStruc(ISAP::PatchExtractor3D* _p, ExtractorType _m, size_t npatch=0);

        DataStruc(){}

        /**
          * Return type of manifold for a given channel
           * param[in] inputData: ptr to array
           * param[in] ChannelDim: vector describing each patch (size=nb channels)
           * param[in] ManifoldType: vector describing Manifolds
           * param[in] isContiguous: flag to specify if data on the same manifold are
           contiguous for all patches (@see DataStructure::m_contiguousChannels)
        */
        DataStruc(dblarray* inputData,std::vector<size_t>& ChannelDim,
            std::vector<Manifold>&  ManifoldType, bool isContiguous, bool OwnPtr=false);

        /**
          * Destructor.
        */
        virtual ~DataStruc(){
            if(m_inputPatches!=nullptr) delete m_inputPatches;}
        /**
          * Iterator
        */
        //class ManifoldIterator;
        /**
          * Signature of class.
        */
        virtual std::string getSignature() const ;

        //Getters
        /**
          * Check if patches have already been extracted
           * @return true if patches extracted.
        */
        bool isExtracted() const { return m_extracted;}

        /**
          * Check if patches have already been extracted
           * @return true if patches extracted.
        */
        bool isExtractorSet() const { if(m_extractor!=nullptr) return true; return false;}

        /**
          * Return number of manifolds in product space
           * @return number of manifolds.
        */
        size_t getNManifold() const {
            return m_manifold.size();
        }
        /**
          * Return type of manifold for a given channel
           * param[in] Channel: index of manifold
           * @return type of selected manifold.
        */
        Manifold getManifold(size_t Channel) const {
            if(Channel<getNManifold()) return m_manifold[Channel];
            else {
                std::cout<<"Incorrect Channel"<<Channel<<" vs "<<m_manifold.size()<<
                                                                                 std::endl;
                MyExit::exit(EXIT_FAILURE);
            }
            return Manifold::Rn;
        }
        /**
          * Return Dimension of a given channel
           * param[in] Channel: index of manifold
           * @return dimension of selected manifold.
        */
        size_t getManifoldDimSize(size_t Channel) const {
            if(Channel<getNManifold()) return m_productDim[Channel];
            else {
                std::cout<<"Incorrect Channel"<<Channel<<" vs "<<m_manifold.size()<<
                                                                                 std::endl;
                MyExit::exit(EXIT_FAILURE);
            }
            return 0;
        }
        /**
          * Return Offset in a line for a given channel (useful to locate a sample)
           * param[in] Channel: index of manifold
           * @return offset of selected manifold.
        */
        size_t getManifoldOffset(size_t Channel) const {
            if(Channel<getNManifold()) return m_offsetDim[Channel];
            else {
                std::cout<<"Incorrect Channel"<<Channel<<" vs "<<m_manifold.size()<<
                                                                                 std::endl;
                MyExit::exit(EXIT_FAILURE);
            }
            return 0;
        }

        /**
          * Create an array containing all data for the selected channel
           * param[in] Channel: index of manifold
           * @return array containing all data for the selected manifold.
        */
        dblarray* getManifoldData(size_t Channel);

        /**
          * Return if multichannel data
           * @return multichannel flag.
        */
        bool getMultiChannelFlag() const {
            if(getNManifold()>1) return true;
            else return false;
        }

        /**
          * Return total number of patches extracted
           * @return number of patches extracted.
        */
        size_t getNPatches() const { return m_nPatches;}

        /**
          * Return size of Patch (=potentially 3D)
           * @return multichannel flag.
        */
        size_t getPatchSize() const { return m_PatchSize;}

        /**
          * Return 2D size
           * @return multichannel flag.
        */
        size_t get2DPatchSize() const { return m_2DPatchSize;}

        /**
          * Return total number of channels per patch
           * @return number of channels.
        */
        size_t getpatchSizeZ() const { return m_ZPatchSize;}

        /**
          * Return a pointer to an array containing the patches
           * @return This extracts an 1D array  containing all channels and patches
            (@see DataStruc::m_inputPatches)
        */
        dblarray* getDataPtr(){ return m_inputPatches;}
        /**
          * Return a pointer to an array containing the patch required
           * param[in] sample: index of the patch among set
           * param[in] channel: index of the channel required among the pacth
           * @return ptr to patch required
        */
        inline double* getSampleChannel(size_t sample,size_t channel){
            if(m_contiguousChannels)
                return m_inputPatches->buffer()+sample*m_productDim[channel]
                                                    +m_offsetDim[channel]*m_nPatches;
            else
                return m_inputPatches->buffer()+m_offsetDim[channel]+getPatchSize()*sample;
        }


        //Setters
        /**
          * Reorder the buffer from/to contiguous
           * param[in] ContiguousFlag: flag
           * @brief Switch from channel contiguous order per patch to overall patches
           channel contiguous and vice-versa.
        */
        void setContiguousChannels(bool ContiguousFlag);

        /**
          * Set Description of data (should be before size)
           * @param[in] _manifold: vector containing manifolds.
        */
        void setManifold(std::vector<Manifold> _manifold){
            if(!m_manifold.empty()) m_manifold.clear();
            m_manifold=_manifold;
        };

        /**
          * Set Description of data size
           * @param[in] _dim: dimension of each manifold.
        */
        void setManifoldDim(std::vector<size_t> _dim){
            if(!m_productDim.empty()) m_productDim.clear();
            if(_dim.size()==m_manifold.size()){
                m_productDim=_dim;
            }
            if(!m_offsetDim.empty()) m_offsetDim.clear();
            m_offsetDim.push_back(0);
            for(size_t kprod=1;kprod<_dim.size();++kprod){
                m_offsetDim.push_back(m_offsetDim[kprod-1]+_dim[kprod-1]);
            }
        };

        /**
          * Reset Description of data (should be consistent)
           * @param[in] _manifold: vector containing manifolds.
           * @param[in] _dim: dimension of each manifold.
        */
         void resetManifold(std::vector<Manifold> _manifold, std::vector<size_t> _dim);

        /**
          * Set Number of patches
           * @param[in] _npatches: number of patches to extract/extracted.
        */
        void setNPatches(size_t _npatches) {  m_nPatches=_npatches;}

        /**
          * Set Random number generator
           * @param[in] RNG: random number generator.
        */
        virtual void setRNG(gslRNG* RNG) {
            if(m_extractor==nullptr){
                std::cout<<"No extractor set : no RNG to set."<<std::endl;
                MyExit::exit(EXIT_FAILURE);
            }
            m_RNG=RNG;m_extractor->setRNG(m_RNG);
        }
        /**
          * Set Seed of Random number generator
           * @param[in] seed: seed for random number generator.
        */
        virtual void setSeed(size_t seed) {m_RNG->setSeed(seed);m_seed=seed;}
        /**
          * Reset seed of Random number generator
           * @details This reset the seed to the previous seed applied (default 0)
        */
        virtual void resetSeed() {m_RNG->setSeed(m_seed);}
        /**
          * Set Verbosity
           * @param[in] verb: boolean flag for verbosity.
        */
        virtual void setVerbose(bool verb){
            m_verb=verb;
            if(m_extractor!=nullptr) m_extractor->setVerbose(verb);
        }
        /**
          * Use the patch extractor to extract Patches
          * @return ptr to array containing all data.
           * @details This extracts the number of patches required (@see
            DataStruc::m_nPatches)
        */
        dblarray* extractPatches() ;

        /**
          * Use the patch extractor to gather all Patches
          * @return status of extraction.
        */
        int gatherPatches() ;

        /**
          * Replace the channel selected with data
           * param[in] Data: ptr to array containing new data
           * param[in] Channel: index of manifold to replace
           * @return array containing all data for the selected manifold.
        */
        void putManifoldData(dblarray* Data, size_t Channel);


        /**
          * Replace the channel selected with data
           * param[in] Data: ptr to array containing new data
           * param[in] Channel: index of manifold to replace
           * @return array containing all data for the selected manifold.
        */
        DataStruc* extractContiguousPatches(size_t PatchOffset,size_t NPatches){
            if(NPatches>=m_nPatches) return this;//nothing 
	    if(PatchOffset+NPatches>m_nPatches) NPatches=m_nPatches-PatchOffset;
            DataStruc* newStruc;
            if(!m_contiguousChannels){
                //Had to copy channel by channel, as they are disjoint in memory
                //for same patch
                dblarray* shiftedData=new dblarray(m_PatchSize*NPatches);
                double* ShiftedBuf=shiftedData->buffer();
                for(size_t km=0;km<getNManifold();++km){
                    double* StartIt=getSampleChannel(PatchOffset,km);
                    for(size_t kel=0;kel<NPatches*getManifoldDimSize(km);++kel,
                        ++ShiftedBuf,++StartIt)
                        *ShiftedBuf=*StartIt;
                }
                newStruc= new DataStruc(shiftedData,m_productDim,
                                m_manifold,m_contiguousChannels, true);
            } else {
                dblarray* shiftedData=new dblarray();
                double* StartPtr=m_inputPatches->buffer()+PatchOffset*m_PatchSize;
		std::cout<<"PP="<<PatchOffset*m_PatchSize<<std::endl;
	        std::cout<<"STARTPTR="<<StartPtr<<std::endl;
                shiftedData->alloc(StartPtr,m_PatchSize*NPatches,0,0);
                newStruc= new DataStruc(shiftedData,m_productDim,
                                m_manifold,m_contiguousChannels, true);
            }
	    return newStruc;
        }

};

/*
Alternative would be a pair of iterators I guess
class DataStruc::ManifoldIterator: public std::iterator<std::bidirectional_iterator_tag, double>{
  double* m_ptr;
public:
  MyIterator(double* _x) :m_ptr(_x) {}
  MyIterator(const MyIterator& mit) : m_ptr(mit.m_ptr) {}
  MyIterator& operator++() {p+=m_totalDim;return *this;}
  MyIterator operator++(int) {MyIterator tmp(*this); operator++(); return tmp;}
  bool operator==(const MyIterator& rhs) const {return m_ptr==rhs.m_ptr;}
  bool operator!=(const MyIterator& rhs) const {return m_ptr!=rhs.m_ptr;}
  double& operator*() {return *m_ptr;}
};*/


#endif //DataStructure_H

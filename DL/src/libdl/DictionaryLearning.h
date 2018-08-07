/******************************************************************************
 **                   Copyright (C) 2016 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Florent Sureau [MOD rewritten from S. Beckouche Code]
 **
 **    Date:  06-07/2016
 **
 **    File:  DictionaryLeaning.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces classes for Dictionary learning
 **    -----------
 **  - 06-07/2016: MOD, kSVD implemented
 ******************************************************************************/


#ifndef DictionaryLearning_H
#define DictionaryLearning_H

#include "SparseCoding.h"
#include "OMPCoding.h"
#include "gslRNG.h"
#include "PatchExtractor.hpp"
#include "DictionaryInit.h"
#include "ISAP_verbose.hpp"
#include "myexit.hpp"
#include <list>
/**
 * @class SparseDL
 * @brief Abstract class essentially implementing the basic setters/getters for
   sparse DL.
*/
class SparseDL {
protected:
    SparseCoding *m_coder;                  //!< Ptr to sparse coder, not owned
    DictionaryInit *m_dictionaryInit = NULL;//!<Initial Dictionary
    DataStruc *m_training = NULL;      //!< Ptr to Training set, not owned
    dblarray *m_dictionary;                 //!< Ptr to dictionary, owned
    dblarray *m_codetraining;               //!< Ptr to training set codebook, owned
    dblarray *m_approx;                     //!< Ptr to training sparse approximation, owned
    dblarray *m_l2cost;                     //!< Ptr to cost function, owned
    dblarray *m_patchMean;                  //!< Ptr to mean of each patch, owned
    double m_eps;                           //!< Accuracy parameter
    bool m_verbose;                         //!< Verbosity flag
    size_t m_natoms;                        //!< Number of atoms
    size_t m_npix;                          //!< Number of pixels per atom/training patch
    size_t m_ntrain;                        //!< Number of examples in the training set
    bool m_patchCentering = false;          //!< If true, all patches are centered by mean subtraction
    bool m_trainingSetCentering = false;    //!< If true, the training set is globally centered by subtraction of the mean sample
    unsigned int m_atom_thresh;             //!< minimal number of time an atom should be used
    unsigned int m_niter;                   //!< Number of iterations in the Training
    gslRNG* m_RNG;                          //!< Random number generator, not owned
    size_t m_channel;/**< Index of channel to process in input data.*/
public:
    //Constructors/Destructors
    SparseDL(){init();}
    /**
      * Constructor of SparseDL setting the coder.
       * @param[in] Coder for sparse DL.
     */
    SparseDL(SparseCoding &Coder){init();setCoder(Coder);}
    virtual ~SparseDL(); //!< Default destructor
    void init();

    //Setters
    /**
      * Set Coder.
       * @param[in] coder for sparse DL.
     */
     virtual void setCoder(SparseCoding &coder){
         m_coder=&coder;
         std::cout<<"BASE CODER:"<<m_coder<<std::endl;
         m_coder->setChannel(0);
     }
     /**
      * Set Training set.
       * @param[in] training Training set.
       * @param[in] Channel channel to be processed.
       * @param[in] online Flag to indicate if patches should be extracted online
       * @return flag to indicate if successfull
     */
      virtual int setTraining(DataStruc &training, int Channel=0,
                                                            bool online=false);
     /**
      * Set Training set for online patch extraction.
       * @param[in] patchExtractor Ptr to 3D patch extractor.
       * @param[in] nTrain Number of training samples desired.
       * @return flag to indicate if successfull
     */
     //virtual int setTraining(ISAP::PatchExtractor* patchExtractor,
     //                                                           size_t nTrain);
      /**
      * Set Verbosity.
       * @param[in] verbose flag for verbosity.
     */
      virtual void setVerbose(bool verbose) {m_verbose =verbose;}
     /**
      * Set Number of interations for DL.
       * @param[in] Niter number of iterations.
     */
      virtual void setNIter(unsigned int Niter){m_niter=Niter;}
     /**
      * Set accuracy for epsilon comparison.
       * @param[in] eps accuracy - for float comparison for instance.
     */
      virtual void setEps(double eps){m_eps=eps;}
     /**
      * Set threshold to identify active atoms.
       * @param[in] atom_thresh min number of examples so that atom is active.
     */
      virtual void setAtomThresh(unsigned int atom_thresh){
                                                    m_atom_thresh=atom_thresh;}
     /**
      * Set Random number generator.
       * @param[in] RNG random number generator.
     */
      virtual void setRNG(gslRNG* RNG) {m_RNG=RNG;}

    /**
     * Set Dictionary.
      * @param[in] Dict  Dictionary to set.
    */
    virtual void setDictionary(dblarray const &Dict);


    /**
     * Set Initial Dictionary according to type.
      * @param[in] nAtoms Number of atoms in dcictionary.
      * @param[in] density density for random initialization.
      * @param[in] InitDictionaryType type of dictionary.
    */
    virtual int initDictionary(size_t nAtoms, double density,
                                             InitDictionaryType dictionaryType);
    /**
     * Indicates patch mean will be set to 0.
      * @param[in] patchCentering flag indicating centering of patches.
    */
    virtual void setPatchCentering(bool patchCentering = true);
    /**
     * Indicates mean patch will be set to 0.
      * @param[in] trainingSetCentering
    */
    virtual void setTrainingSetCentering(bool trainingSetCentering = true);

    /**
     * Check if patch mean will be set to 0.
      * @return patchCentering flag indicating centering of patches.
    */
    bool getPatchCentering()const {return m_patchCentering;}
    /**
     * Check if mean patch will be set to 0.
      * @return trainingSetCentering
    */
    bool getTrainingSetCentering()const {return m_trainingSetCentering;}


    //Getters
    /**
     * Get number of Iterations for DL.
      * @return number of DL iterations.
    */
    unsigned int getNIter() const{return m_niter;}
    /**
     * Get DL Dictionary.
      * @return Ptr to Dictionary.
    */
    dblarray* getDictionary() const{return m_dictionary;}
     /**
      * Get Training set.
       * @return Ptr to Training set.
     */
    DataStruc* getTraining() const{return m_training;}

    /**
     * Get Mean of each Patch.
      * @return Ptr (maybe NULL) to array containing the mean of each patch.
    */
    dblarray* getPatchMean() const {return m_patchMean;}

    /**
     * Get Sparse Coder.
      * @return Ptr to Sparse Coder.
    */
    SparseCoding* getCoder() const{return m_coder;}
    /**
     * Get Number of atoms in the Dictionary.
      * @return Number of atoms.
    */
    size_t getNAtoms() const {return m_natoms;}
    /**
     * Get number of pixels in atoms/input training or test set.
      * @return number of pixels.
    */
    size_t getNPix() const {return m_npix;}
    /**
     * Get number of training examples.
      * @return number of training examples.
    */
    size_t getNTrain() const {return m_ntrain;}
    /**
     * Get Signature of initial dictionary
      * @brief SparseDL::getInitialDictionaryName
      * @return Initial dictionary signature
    */
    std::string getInitialDictionaryName() const{
      return m_dictionaryInit->getSignature();}
    /**
     * Interface to signature.
      * @return signature of class.
    */
    virtual std::string getSignature()      const = 0;

    //Processing functions
    /**
      * Interface to Dictionary update.
     */
    virtual void updateDictionary()=0;
    /**
      * Interface to Dictionary learning.
     */
    virtual int  learning()=0;
    /**
      * Interface to dictionary centering (mean of atoms=0).
     */
    virtual void centerDictionary()=0;
    /**
      * Interface to patch centering (mean of patches=0).
     */
    virtual void centerPatches()=0;
    /**
      * Interface to training centering (mean patch over training patches=0).
     */
    virtual void centerTraining(char *nameMean=NULL)=0;
    /**
      * Interface to l2 ball projection of atoms.
     */
    virtual void l2projDictionary()=0;
    /**
      * Interface to l2 ball projection of training patches.
     */
    virtual void l2projTraining()=0;
    /**
      * Interface to centering and l2 ball projection of atoms.
     */
    virtual void normDictionary()=0;
    /**
      * Interface to centering and l2 ball projection of atoms.
     */
    virtual void normTraining()=0;
private:
};

/**
 * @class SparseVectorDL
 * @brief Abstract class essentially implementing the matrix operations and
 getters/setters to perform vector(ized) sparse DL.
*/
class SparseVectorDL: public SparseDL {
protected:
    MatrixDecomposition *m_matdecomp;/**< Matrix routines, not owned.*/
    Manifold m_manifold;
public:
    //Constructors/Destructors
    /**
      * Default Constructor of SparseVectorDL .
     */
    SparseVectorDL(){init();}
    /**
      * Constructor of SparseVectorDL setting the coder.
       * @param[in] Coder for sparse DL.
     */
    SparseVectorDL(SparseCoding &Coder):SparseDL(Coder){init();setCoder(Coder);}
    /**
      * Default Destructor of DLParameters.
     */
    virtual ~ SparseVectorDL(){ }
    /**
      * Initialization of matrix decomposition.
     */
    void init(){m_matdecomp=NULL;m_manifold=Manifold::Rn;}

    /**
      * Initialization of learning.
     */
    int initDL();
    //Getters
    /**
      * Get Training set Code.
       * @param[in] recompute Flag to indicate if sparse coding is performed or
                    if the last code computed is desired.
       * @return Ptr to the sparse code of the training set
     */
    inline dblarray*& computeCode(bool recompute=true){
        if((recompute)||(m_codetraining==NULL)){
            if(m_codetraining!=NULL) delete m_codetraining;
            m_coder->setDictionary(*m_dictionary);
            m_codetraining= m_coder->computeCode(m_training);
        }
        return m_codetraining;
    }

    /**
      * Get Matrix Decomposition associated with learner.
       * @return Ptr to matrix decomposition
     */
    MatrixDecomposition* getMatDecomp() const { return m_matdecomp;}


    /**
        * Get sparse approximation.
         * @param[in] recompute Flag to indicate if sparse coding is performed or
                    if based on the last code computed.
         * @param[in] overall Flag to indicate if mean added back in case of
                      patch centering.
         * @return Ptr to the sparse approximation of the training set
    */
    inline dblarray*& getApprox(bool recompute=true, bool overall=false) {
        if((recompute)||((m_codetraining==NULL)&&(m_approx==NULL))){
            m_codetraining = computeCode(recompute);
            if(m_approx != NULL) delete m_approx;
            m_approx = NULL;
        }
        if(m_approx == NULL){
            //m_coder->setDictionary(*m_dictionary);
            m_approx=m_coder->computeApprox(m_codetraining,m_training);
        } //m_approx=m_matdecomp->matMultTransp(*m_codetraining,*m_dictionary);
        if((overall) && (m_patchCentering)){//We want to add the mean back
            if(m_coder->getManifold()==Manifold::Rn)
                m_matdecomp->addAlongDim(*m_approx, 0, *m_patchMean,1.0);
            else {
                std::cout<<"APPROX: DO NOT ADD MEAN BACK"<<std::endl;
            }
        }
        return m_approx;
    }

    /**
      * Get Manifold.
       * @return Manifold where DL operates on
     */
    Manifold getManifold() const {return m_manifold;};

    //Setters

    /**
      * Set Coder.
       * @param[in] coder for sparse DL.
     */
     virtual void setCoder(SparseCoding &coder){
         std::cout<<"DERIVED CODER:"<<m_coder<<std::endl;
         setManifold(coder.getManifold());
     }

    /**
      * Set Channel to process for input data structure.
       * @param[in] dictionary dictionary to set.
    */
    virtual void setChannel(const size_t _chan=0){
        m_channel=_chan;
        m_coder->setChannel(_chan);
    }
    /**
      * Set Manifold DL operates on.
       * @param[in] _m Manifold.
    */
    void setManifold(Manifold _m) { m_manifold=_m;};


    /**
      * Set initial dictionary for learning.
       * @param[in] initDict initial dictionary.
     */
    virtual void setDictionary(dblarray &initDict); /// @todo AW readd const
    /**
     * Set Training set.
      * @param[in] training Training set.
      * @param[in] Channel channel to be processed.
      * @param[in] online Flag to indicate if patches should be extracted online
      * @return flag to indicate if successfull
    */
     virtual int setTraining(DataStruc &training, int Channel=0,
                                                           bool online=false);
    /**
      * Set matrix operation classes.
       * @param[in] mdecomp matrix operation class.
     */
    void setMatDecomp(MatrixDecomposition &mdecomp) {m_matdecomp=&(mdecomp);}
    //Implementation of pre-processing for matrices/vectors
    /**
      * Center the dictionary (mean over pixels in each atom =0).
     */
    virtual void centerDictionary();
    /**
      * Center the patches (mean over pixels in each patch =0).
     */
    virtual void centerPatches();
    /**
      * Center the training set (mean over patches in each pixel =0).
     */
    virtual void centerTraining(char *nameMean=NULL);
    /**
      * Project each atom on the l2 ball.
     */
    virtual void l2projDictionary();
    /**
      * Project each patch on the l2 ball.
     */
    virtual void l2projTraining();
    /**
      * Center and project each atom on the l2 ball.
     */
    virtual void normDictionary();
    /**
      * Center and project each atom on the l2 ball.
     */
    virtual void normTraining();
    /**
      * Replace atoms not sufficiently used (@see SparseDL::m_atom_thresh).
     */
    unsigned int replaceUnusedAtoms();
    /**
      * Perform learning.
     */
    virtual int learning();
    /**
      * Replace non significant or highly correlated (=) atoms
     */
    unsigned int cleanDictionary();
    /**
      * Compute average sparse approximation error
       * @param[in] recompute flag to recompute the code.
     */
    double computeAverageError(bool recompute);

protected:
    /**
      * Get given atom in the dictionary
       * param[in] AtomNb index of the atom in the columns of the dictionary
       * param[out] Atom Ptr to atom values
     */
    inline void getAtom(size_t AtomNb,double *Atom){
        double *iterBuf= m_dictionary->buffer()+ AtomNb;
        double *iterAtom=Atom;
        for(size_t np=0;np<m_npix;++np, iterBuf +=m_natoms,++iterAtom)
            *iterAtom=*iterBuf;
    }
    /**
      * Get given example in the training set
       * param[in] sampleNb index of the example in the training set
       * param[out] TrainSample Ptr to patch values
     */
    inline void getTrainingSample(size_t sampleNb,double *TrainSample){
        m_coder->getInputSample(m_training,sampleNb,TrainSample);
    }
    /**
      * Get given example in the training set
       * param[in] sampleNb index of the example in the training set
       * param[out] TrainSample Ptr to patch values
     */
    inline void getTrainingCode(size_t sampleNb,double *SampleCode){
        double *iterBuf= m_codetraining->buffer()+ sampleNb* m_natoms;
        double *iterCode=SampleCode;
        for(size_t np=0;np<m_natoms;++np, ++iterBuf,++ iterCode)
            *iterCode =*iterBuf;
    }
    /**
      * Set an atom in the dictionary
       * param[in] AtomNb index of the atom in the dictionary (column)
       * param[out] Atom Ptr to atom values
       * param[in] stride to get consecutive pixels in the atom
     */
    inline void setAtom(size_t AtomNb,double *Atom,size_t stride){
        double *iterBuf= m_dictionary->buffer()+ AtomNb;
        double *iterAtom=Atom;
        for(size_t np=0;np<m_npix;++np, iterBuf +=m_natoms,iterAtom+=stride)
            *iterBuf=*iterAtom;
    }
    /**
      * Set a training sample according to its position
       * param[in] sampleNb index of the example (row).
                   Pixels are contiguous in memory.
       * param[out] TrainSample Ptr to training patch values
     */
    inline void setTrainingSample(size_t sampleNb,double *TrainSample){
        double *iterBuf= m_training->getSampleChannel(sampleNb,m_channel);
        double *iterTrain= TrainSample;
        for(size_t np=0;np<m_npix;++np, ++iterBuf,++iterTrain)
            *iterBuf =*iterTrain;
    }
    /**
      * Set a training code according to its position
       * param[in] sampleNb index of the example (row).
                   Pixels are contiguous in memory.
       * param[out] SampleCode Ptr to the sparse code of the example.
     */
    inline void setTrainingCode(size_t sampleNb,double *SampleCode){
        double *iterBuf= m_codetraining->buffer()+ sampleNb* m_natoms;
        double *iterCode=SampleCode;
        for(size_t np=0;np<m_natoms;++np, ++iterBuf,++ iterCode)
            *iterBuf =*iterCode;
    }

};

/**
 * @class SparseVectorMOD
 * @brief Class implementing the method of optimal directions (MOD), an
 * alternated minimization scheme on code and dictionary update where the
 * dictionary update is performed in one step using the code pseudo-inverse.
*/
class SparseVectorMOD: public SparseVectorDL {
public:
    /**
      * Default constructor of SparseVectorMOD.
     */
    SparseVectorMOD(): SparseVectorDL() {m_grad_descent=false;}
    /**
      * Constructor of SparseVectorMOD setting the coder.
       * @param[in] Coder for sparse DL.
     */
    SparseVectorMOD(SparseCoding &Coder, bool _grad=false):
                SparseVectorDL(Coder){}


    /**
      * Default Destructor of SparseVectorMOD.
     */
    virtual ~SparseVectorMOD(){ }
    /**
      * Implementation of the signature
       *@return Signature of SparseVectorMOD
     */
    virtual std::string getSignature() const {
        return std::string("Sparse Vectorized MOD");}

    /**
      * Implementation of Dictionary Update once sparse coding is performed
     */
    virtual void updateDictionary();
protected:
    bool m_grad_descent;//!< Flag to indicate gradient descent is performed
    size_t m_niter_grad;//!< Max number of gradient descent iterations
    size_t m_nlsmax;//!< Max number of iterations for step selection
    double m_stopGrad;//!<  Stopping criterion based on max abs gradient value
};



/**
 * @class SparseVectorKSVD
 * @brief Class implementing two K-SVD variants, an alternated minimization
  * scheme where after sparse coding, a dictionary is learned atom by atom by
  * looking at the best rank-1 approximation of the partial residuals for
  * examples using this atom [subtracting the contributions of all other atoms].
  * This can be achieved either using SVD or an iterative power method.
*/
class SparseVectorKSVD: public SparseVectorDL {
protected:
    bool m_approxSVD;/**< Flag to specify if approximate K-SVD or K-SVD.*/
    unsigned int m_approxSVDIter;/**< Number of iterations in the power method.*/
public:
    /**
      * Default Constructor of SparseVectorKSVD.
     */
    SparseVectorKSVD(): SparseVectorDL(), m_approxSVD(false),
        m_approxSVDIter(1){}
    /**
      * Constructor of SparseVectorKSVD setting the coder.
       * @param[in] Coder for sparse DL.
      * @brief by default, this set SparseVectorKSVD::m_approxSVD to false and
      * SparseVectorKSVD::m_approxSVDIter to 1.
     */
    SparseVectorKSVD(SparseCoding &Coder): SparseVectorDL(Coder),
        m_approxSVD(false),
        m_approxSVDIter(1){}
    /**
      * Constructor of SparseVectorKSVD setting the coder, K-SVD or approximate
        variant.
       * @param[in] Coder for sparse DL.
       * @param[in] approx Flag to set approximate K-SVD.
       * @param[in] iter Number of iterations in power method if approximate
       *            K-SVD is performed (default 1).
     */
    SparseVectorKSVD(SparseCoding &Coder,bool approx,unsigned int iter=1):
        SparseVectorDL(Coder){
        setApproxKSVD(approx,iter);
    }
    /**
      * Default Destructor of SparseVectorKSVD.
     */
    virtual ~ SparseVectorKSVD(){ }
    /**
      * Set which K-SVD variant used.
       * @param[in] approx Flag to set approximate K-SVD.
       * @param[in] iter Number of iterations in power method if approximate
       *            K-SVD is performed (default 1).
     */
    void setApproxKSVD(bool approx,unsigned int iter=1){
        m_approxSVD =approx;
        m_approxSVDIter =iter;
    }
    /**
      * Implementation of signature
       * @return Signature of SparseVectorKSVD.
     */
    virtual std::string getSignature() const {
        if(m_approxSVD) return std::string("Sparse Vectorized K-SVD");
        else return std::string("Sparse Vectorized approximate K-SVD");
    }
    /**
      * Implementation of dictionary update step for K-SVD
     */
    virtual void updateDictionary();
protected:
    /**
      * Create a list of samples using a given atom
       * @param[in] AtomNb Index of atom to consider
       * @param[out] listTrain Ptr to list of samples using the atom.
       * @param[out] listCoefs Ptr to values of the atom coefs for samples in
                     the list.
       * @return Number of elements in the list
     */
    size_t createListSubMatrix(size_t AtomNb,size_t*listTrain,
                               double* listCoefs);
};


#endif //DictionaryLearning_H

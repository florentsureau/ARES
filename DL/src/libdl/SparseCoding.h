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
 **    File:  SparseCoding.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces abstract classes for Sparse decomposition
 **    -----------
 ** - 06/2016: rewriting code with abstract classes
 **
 ******************************************************************************/

#ifndef SparseCoding_H
#define SparseCoding_H

#include "DataStructure.h"

/**
 * @class SparseCoding_Settings
 * @brief Settings for SparseCoding class.
*/
class SparseCoding_Settings {
public:
    bool verbose;/**< Verbosity flag.*/
    double errorTarget;/**< L2-norm approximation error targeted.*/
    size_t sparsityTarget;/**< Max number of atoms active for each example.*/
    size_t minimalSparsity;/**< Min number of atoms active for each example.*/
    double minimalCorrelation; //!< Minimum correlation value to validate an atom selection (for one pixel)
    double meanSparsityTarget;//!< target average sparsity
    double epsilon;/**< Precision used for comparisons.*/
    Manifold m_manifold;
    SparseCoding_Settings();
    ~SparseCoding_Settings() {}

    void init();
    void print(string sparseCodingName = "") const;
};


/**
 * @class SparseCoding
 * @brief Abstract class for sparse coding.
*/
class SparseCoding {
protected:
    SparseCoding_Settings m_settings;/**< Settings for Sparse Coding.*/
    dblarray *m_dictionary;     /**< Ptr to dictionary (nAtoms*nPix), not owned.*/
    size_t m_natoms;            /**< Number of atoms in dictionary.*/
    size_t m_npix;              /**< Number of samples per atom.*/
    size_t m_ninput;            /**< Number of vectors to sparse code.*/
    dblarray m_rfactor;         /**< rescaling factor for code selection.*/
    dblarray m_inputMetric;     /**< Metric to use in input space.*/
    dblarray m_codeMetric;      /**< Metric to use in code space.*/
    bool m_boolInputMetric;     /**< Flag set if input metric for atom selection.*/
    bool m_boolDiagInputMetric; /**< Flag set if input metric diagonal.*/
    bool m_boolCodeMetric;      /**< Flag set if code metric is used. //not implemented yet*/
    bool m_boolDiagCodeMetric;  /**< Flag set if code metric diagonal.*/

    bool m_matdecompIsLocal;/**< Matrix decomposition owned if true.*/
    MatrixDecomposition *m_matdecomp = nullptr;/**< Ptr to matrix operations.*/
    bool m_isAllocated = false;/**< Flag for Dictionary/decomposition allocation.*/
    /**
      * Check allocation status of Sparse Coder.
      * @return EXIT_SUCCESS, EXIT_FAILURE depending on allocation status
     */
    virtual int checkAllocation();

    size_t m_channel=0;/**< Index of channel to process in input data.*/

public:
    //constructors/destructors
    /**
      * Standart constructor.
     */
    SparseCoding();
    /**
      * Constructor with dictionary included.
       * @param[in] dictionary set dictionary used for sparse coding
     */
    SparseCoding(dblarray &dictionary);
    /**
      * Standart destructor.
    */
    virtual ~SparseCoding();
    /**
      * Default initialization of members in the class.
    */
    void init();

    //i/o
    /**
      * Read input metric from fits file.
       * @param[in] FitsFileName filepath to fits with input metric.
    */
    dblarray* readMetric(const char *FitsFileName);
    /**
      * Print settings for sparse coding @see SparseCoding_Settings
    */
    void print();
    virtual std::string getName()=0;
    //setters
    /**
      * Set dictionary.
       * @param[in] dictionary dictionary to set.
    */
    virtual void setDictionary(const dblarray &dictionary);
    /**
      * Set Channel to process for input data structure.
       * @param[in] dictionary dictionary to set.
    */
    virtual void setChannel(const size_t _chan){m_channel=_chan;}

    /**
      * Set code rescale factor (used for selection of atoms).
       * @param[in] dictionary dictionary to set.
    */
    virtual void setRFactor(const dblarray &rescale_factor);
    /**
      * Interface to set input metric.
       * @param[in] metric input metric to set.
    */
    virtual void setInputMetric(const dblarray &metric)=0;
    /**
      * Interface to set code metric.
       * @param[in] metric code metric to set.
    */
    virtual void setCodeMetric(const dblarray &metric)=0;
    /**
      * Interface to set matrix decompositions.
       * @param[in] mdecomp matrix decomposition to use.
    */
    virtual void setMatDecomp(MatrixDecomposition &mdecomp)=0;
    /**
      * Interface to set matrix decompositions.
    */
    virtual void setMatDecomp()=0;
    /*!
     * \brief setParameters
     * \param i_setting
     */
    virtual void setParameters(const SparseCoding_Settings &i_setting) {
        m_settings = i_setting;
    }
    /**
      * Set targeted l2-norm of approximation residuals.
       * @param[in] errortarget l2-norm targeted.
    */
    virtual void setErrorTarget(const double errortarget) {
      m_settings.errorTarget= errortarget;}
    /**
      * Set targeted sparsity level of code.
       * @param[in] sparsitytarget sparsity level of code.
    */
    virtual void setSparsityTarget(const size_t sparsitytarget) {
      m_settings.sparsityTarget= sparsitytarget;}
    /**
      * Set targeted minimal sparsity level of code.
       * @param[in] minimalSparsity minimal sparsity level of code.
    */
    virtual void setMinimalSparsity(const size_t minimalSparsity) {
      m_settings.minimalSparsity= minimalSparsity;}
    /**
      * Set targeted minimal correlation level of atom.
       * @param[in] minimalCorrelation minimal correlation level of atom.
    */
    virtual void setMinimalCorrelation(const double minimalCorrelation) {
      m_settings.minimalCorrelation = minimalCorrelation;}
    /**
      * Set verbosity.
       * @param[in] verbose flag for verbosity.
    */
    void setVerbose(bool verbose){
      m_settings.verbose = verbose;}
    /**
      * Set precision parameter.
       * @param[in] eps precision parameter.
    */
    void setEps(double eps){m_settings.epsilon=eps;}
    //getters
    /**
      * Get ptr to dictionary.
       * @return ptr to dictionary.
    */
    virtual dblarray* getDictionary() const {return m_dictionary;}
    /**
      * Get number of atoms in dictionary.
       * @return number of atoms in dictionary.
    */
    virtual Manifold getManifold() const { return m_settings.m_manifold;}
    /**
      * Get number of atoms in dictionary.
       * @return number of atoms in dictionary.
    */
    size_t getNbAtoms() const { return m_natoms;}
    /**
      * Get number of samples per atom.
       * @return number of samples per atom.
    */
    size_t getNbpixels() const { return m_npix;}
    /**
      * Get sparsity level targeted.
       * @return max number of active atoms per sample.
    */
    size_t getSparseTarget() const { return m_settings.sparsityTarget;}
    /**
      * Get approximation error targeted.
       * @return max l2 norm of approximation residuals.
    */
    double getErrorTarget() const { return m_settings.errorTarget;}

    /**
      * Get an input sample from its index .
       * @param[in] InData input Data
       * @param[in] sampleNb index of the sample
       * @param[out] InSample ptr to buffer where the result is stored.
       * @param[in] Channel index of channel to get
    */
    virtual void getInputSample(DataStruc*& InData,size_t sampleNb,
                                                    double *InSample) =0;

    /**
      * Interface to get the reference tangent plane if it exists.
       * @return ptr to array with reference tangent plane, or nullptr if it is
       not relevant.
    */
    virtual dblarray* getRefTangentPlane() const {
        std::cout<<"NO REF TANGENT PLANE"<<std::endl;
        return nullptr;
        };

    /**
      * Interface to get the signature of the object.
       * @return signature of the object.
    */
    virtual std::string getSignature()=0;
    /**
      * Return the params for sparse coding.
       * @return SparseCoding_Settings object.
    */
    virtual SparseCoding_Settings getParam() const{return m_settings;}

    //processing functions
    /**
      * Apply square root of weights to the data and dictionary.
       * @param[in] Input data to weight.
       * @param[in] Weights metric used for WLS.
       * @return input weighted by square root of weights
       * @note the dictionary is also weighted accordingly
    */

    virtual dblarray* applySqrtWeights(dblarray*& Input, dblarray& Weights)=0;

    /**
      * Interface to get the code of the samples.
      * @param[in] InputData: ptr to input structure.
       * @return code of the samples.
    */
    virtual dblarray* computeCode(DataStruc* InputData)=0;

    /**
      * Interface to get the sparse approximation from input codes.
       * @param[in] Code ptr to input code.
       * @param[in] Input ptr to input structure.
       * @return ptr to array with sparse approximation.
    */
    virtual dblarray* computeApprox(dblarray* Code ,DataStruc* Input)=0;

    /**
      * Interface to get the sparse approximation from input codes.
       * @param[in] Code ptr to input code.
       * @param[out] Approx ptr to output sparse approximation.
       * @param[in] Input ptr to input structure.
       * @return ptr to array with sparse approximation.
    */
    virtual void computeApprox(dblarray* Code, dblarray* Approx, DataStruc* Input)=0;




protected:
    /**
      * Get a pointer on the dictionary buffer.
       * @return pointer on the dictionary buffer.
    */
    double *getDictionaryBuffer(){ return m_dictionary->buffer();}
};


#endif //SparseCoding_H

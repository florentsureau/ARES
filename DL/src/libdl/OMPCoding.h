/******************************************************************************
 **                   Copyright (C) 2016 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Florent Sureau [Basic OMP rewritten from S. Beckouche Code]
 **
 **    Date:  06-07/2016
 **
 **    File:  OMPCoding.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces OMP variants  for Sparse decomposition
 **    -----------
 ** - 06/2016: OMP added
 ** - 01-02/2017: OMP for S1 data added
 ** - 06/2017: Product Space Sparse coding
 ******************************************************************************/

#ifndef OMPCODING_H
#define OMPCODING_H

#include "SparseCoding.h"
#include "OpenMPInterface.h"
#include <limits>

/**
 * @class OMP
 * @brief Class implementing orthogonal matching pursuit (OMP).
*/
class OMP : public SparseCoding {
protected:
    dblarray* m_weightedDictionary;/**< Ptr to dictionary (nAtoms*nPix),
                                weighted by sqrt of input metric,  owned.*/
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    OMP():SparseCoding(),m_weightedDictionary(nullptr){
        m_settings.m_manifold=Manifold::Rn;
    }
    /**
      * Standart destructor.
    */
    virtual ~OMP();


    /**
      * AW addition.
    */
    virtual std::string getName() { return std::string("OMP"); }

    //setters
    /**
      * Set dictionary, overriding base class implementation.
       * @param[in] dictionary dictionary to set.
    */
    virtual void setDictionary(const dblarray &dictionary);
    /**
      * Set code rescale factor, overriding base class implementation.
       * @param[in] dictionary dictionary to set.
    */
    virtual void setRFactor(const dblarray &rescale_factor);
    /**
      * Set input metric, overriding base class implementation.
       * @param[in] metric input metric to set.
    */
    virtual void setInputMetric(const dblarray &metric);
    /**
      * Set code metric, overriding base class implementation.
       * @param[in] metric code metric to set.
    */
    virtual void setCodeMetric(const dblarray &metric);
    /**
      * Set matrix decomposition, overriding base class implementation.
       * @param[in] mdecomp decomposition to use.
    */
    virtual void setMatDecomp(MatrixDecomposition &mdecomp);
    /**
      * Set matrix decomposition, overriding base class implementation.
    */
    virtual void setMatDecomp();

    //getters
    /**
      * Implement signature of OMP.
       * @return signature of OMP
    */
    virtual std::string getSignature();

    /**
      * Get an input sample from its index .
       * @param[in] InData input Data
       * @param[in] sampleNb index of the sample
       * @param[out] InSample ptr to buffer where the result is stored.
    */
    virtual void getInputSample(DataStruc*& InData,size_t sampleNb,
                                            double *InSample);


    //processing functions


    /**
      * Apply square root of weights to the data and dictionary.
       * @param[in] Input data to weight.
       * @param[in] Weights metric used for WLS.
       * @return input weighted by square root of weights
       * @note the dictionary is also weighted accordingly
    */

    dblarray* applySqrtWeights(dblarray*& Input, dblarray& Weights);

    /**
      * Get input space inner product.
       * @param[in] vector1 first input space vector.
       * @param[in] vector2 second input space vector.
       * @return input space inner product.
    */
    virtual double inputDotProduct(const dblarray &vector1,
                                   const dblarray &vector2);
    /**
      * Get input space norm of a vector.
       * @param[in] vector1 input space vector.
       * @return input space norm.
    */
    virtual double inputNorm(const dblarray &vector1);
    /**
      * Get code space inner product.
       * @param[in] vector1 first code vector.
       * @param[in] vector2 second code vector.
       * @return code space inner product.
    */
    virtual double codeDotProduct(const dblarray &vector1,
                                  const dblarray &vector2);
    /**
      * Get correlation in between input space example and an atom
       * @param[in] ve1 first input space vector.
       * @param[in] atom atom of the dictionary.
       * @param[in] inmetric metric used for inner product.
       * @return coorelation between atom and input sample.
    */
    virtual double atomsDotProduct(double *vec1, double *atom, double *inmetric);
    /**
      * Get closest atom to the residual (maximum correlating)
       * @param[in] residualPtr ptr to the residual.
       * @param[in/out] codePtr ptr to array where correlations are stored.
       * @param[in] inmetric metric used for inner product.
       * @param[in] scalePtr additional scaling of the atoms.
       * @return atom index with maximal correlation.
    */
    virtual int findClosestAtom(double* residualPtr,double* codePtr,
                                      double *inputMetricPtr, double* scalePtr);
    /**
      * Get code space norm of a vector.
       * @param[in] vector1 code vector.
       * @return code space norm.
    */
    virtual double codeNorm(const dblarray &vector1);
    /**
      * Implementation to get codes for input examples using OMP.
       * @param[in] InputData ptr to input structure.
       * @return input set code .
    */
    virtual dblarray* computeCode(DataStruc* InputData);

    /**
      * Implementation to get the sparse approximation from input codes.
         * @param[in] Code ptr to input code.
         * @param[in] Input ptr to input structure.
         * @return ptr to array with sparse approximation.
    */
    virtual dblarray* computeApprox(dblarray* Code,DataStruc* Input=nullptr);

    /**
      * Interface to get the sparse approximation from input codes.
       * @param[in] Code ptr to input code.
       * @param[out] Approx ptr to output sparse approximation.
       * @param[in] Input ptr to input structure.
       * @return ptr to array with sparse approximation.
    */
    virtual void computeApprox(dblarray* Code, dblarray* Approx, DataStruc* Input);
};

#endif //OMPCODING_H

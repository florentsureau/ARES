/******************************************************************************
 **                   Copyright (C) 2017 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Florent Sureau
 **
 **    Date:  09/2017
 **
 **    File:  GradOperators.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces Gradient Operator virtual class
 **    -----------
 ** - 09/2017: Virtual Class with gradient properties
 ******************************************************************************/

#ifndef GRADOP_H
#define GRADOP_H

#include "MatrixDecomposition.h"
#include "S1Decomposition.h"
#include "OpenMPInterface.h"

#include <limits>
#include <vector>
#include <list>
#include <algorithm>
#include <functional>



/**
 * @class GradOperator
 * @brief Interface class for gradient operators.
*/
template <typename PARAM>  class GradOperator {
protected:
    PARAM m_precision;/**<Precision.*/
    PARAM m_lipschitz;/**<Lipschitz constant if applicable.*/
    PARAM m_lipschitzFactor;/**< Upper bound factor used to compute Lipschitz value.*/
public:
    //constructors/destructors
    /**
      * Standart constructor
    */
    GradOperator(){init();}

    /**
      * Standart destructor.
    */
    virtual ~GradOperator(){}

    //Setters
    /**
      * Default initialization.
    */
    void init();
    /**
      * Default precision.
    */
    void setDefaultPrecision();


    //Getters
    /**
      * Interface to get the signature of the object.
       * @return signature of the object.
    */
    virtual std::string getSignature() const =0 ;

    /**
      * Get Lipschitz constant for the gradient, if applicable
       * @return Lipschitz constant .
    */
    PARAM getGradLipschitzConstant() const {return m_lipschitz;}

    /**
      * Get Lipschitz constant for the gradient composed with affine transform,
        if applicable, or upper bound.
       * @return Lipschitz constant .
    */
    virtual PARAM getGradLipschitzConstant(const PARAM* LinearTrans,const size_t Nrows,
        const size_t &Ncols) const =0;

    /**
      * Interface to compute the gradient of the function.
       * @param[in] RefData point where to compute the gradient
       * @param[in] Nelems number of elements in reference point
       * @return array containing the gradient.
    */
    virtual PARAM* computeGradient(const PARAM* RefData, const size_t &Nelems)=0;


    /**
      * Interface to compute the gradient of a function composed with a linear
       operator f(Ax+b).
       * @param[in] RefData point where to compute the gradient
       * @param[in] OffsetVal offset point b
       * @param[in] LinearTrans linear transform A
       * @param[in] Nrows number of observations
       * @param[in] Ncols number of elements in reference point
       * @return array containing the gradient.
    */
    PARAM* computeAffineTransformGradient(const PARAM* RefData,
                                const PARAM* OffsetVal,const PARAM* LinearTrans,
                                        const size_t Nrows,const size_t &Ncols);


    /**
      * Compute value of function at given point.
       * @param[in] RefData point where to compute the function
       * @param[in] Nelems number of elements in reference point
       * @return  value of function at the given point.
      @note The default version is using the Moreau identity. If known, this
                                                           can be computed directly.
                                                        */
    virtual PARAM computeFunction(PARAM* RefData,const size_t &Nelems)=0;

    /**
      * Interface to compute the value of a function composed with a linear
       operator f(Ax+b).
       * @param[in] RefData point where to compute the gradient
       * @param[in] OffsetVal offset point b
       * @param[in] LinearTrans linear transform A
       * @param[in] Nrows number of observations
       * @param[in] Ncols number of elements in reference point
       * @return value of the function.
    */
    PARAM computeAffineTransformFunction(const PARAM* RefData,
                                const PARAM* OffsetVal,const PARAM* LinearTrans,
                                        const size_t Nrows,const size_t &Ncols);

};

/**
 * @class GradL2Fidelity
 * @brief Class for (Weighted) l2 norm on residual.
*/
template <typename PARAM>  class GradL2Fidelity:public GradOperator<PARAM>{
protected:
    PARAM* m_refPoint;/**<Ptr to reference - not owned.*/
    PARAM* m_weight;/**<Ptr to Weighting matrix for L2 fidelity term -not owned.*/
    bool m_weightFlag;/**<Flag for weighting matrix.*/
    size_t m_weightNCols;/**<Number of columns for weight matrix.*/
    size_t m_weightNRows;/**<Number of rows for weights- 0 if diagonal.*/
public:
    //constructors/destructors
    /**
      * Standart constructor
    */
    GradL2Fidelity(const PARAM* ref_point,const size_t& Ncols,
                        const PARAM* weights=NULL,const size_t Nrows=0){
        init();
        m_refPoint=ref_point;
        setWeights(weights,Ncols,Nrows);
    }

    /**
      * Standart destructor.
    */
    virtual ~GradL2Fidelity(){}

    //Setters
    /**
      * Default initialization.
    */
    void init();

    /**
      * Set Covariance matrix precision.
    */
    void setWeights(const PARAM* weights,const size_t& Ncols,const size_t& Nrows);

    //Getters
    /**
      * Interface to get the signature of the object.
       * @return signature of the object.
    */
    virtual std::string getSignature() const{
        std::stringstream strstr;
        std::string Sign=std::string("Grad-L2");
        if(m_weightFlag){
            if(m_weightNRows==0) strstr<<"-DiagWeight";
            else strstr<<"-Weight";
            Sign.append(strstr.str());
        }
        return Sign;
    }

    /**
      * Interface to compute the gradient of the function.
       * @param[in] RefData point where to compute the gradient
       * @param[in] Nelems number of elements in reference point
       * @return array containing the gradient.
    */
    virtual PARAM* computeGradient(const PARAM* RefData, const size_t &Nelems);

    /**
      * Compute value of function at given point.
       * @param[in] RefData point where to compute the function
       * @param[in] Nelems number of elements in reference point
       * @return  value of function at the given point.
      @note The default version is using the Moreau identity. If known, this
                                                           can be computed directly.
                                                        */
    virtual PARAM computeFunction(PARAM* RefData,const size_t &Nelems);

    /**
      * Get Lipschitz constant for the gradient composed with affine transform,
        if applicable, or upper bound. Here it is ||A||| |W||
       * @return Lipschitz constant .
    */
    virtual PARAM getGradLipschitzConstant(const PARAM* LinearTrans,const size_t Nrows,
        const size_t &Ncols) const;


};

#endif //GRADOP_H

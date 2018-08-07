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
 **    File:  ProxOperators.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces Proximity Operator virtual class
 **    -----------
 ** - 09/2017: Virtual Class with prox properties
 ******************************************************************************/

#ifndef PROXOP_H
#define PROXOP_H

#include "MatrixDecomposition.h"
#include "S1Decomposition.h"
#include "OpenMPInterface.h"
#include "DataStructure.h"

#include <limits>
#include <vector>
#include <list>
#include <algorithm>
#include <functional>

/**
 * @enum StopCriteria
 * @brief choice among stopping criterion
*/
enum class StopCriteria {Niter,GradLinf,CostCond,RelSol};


/**
 * @class ProxOperator
 * @brief Interface class for proximity operators.
*/
template <typename PARAM>  class ProxOperator {
protected:
    size_t m_nIter;/**<Number of Iterations to compute the prox if iterative.*/
    PARAM m_minGradlinf;/**<Gradient condition if iterative.*/
    PARAM m_minRelCost;/**<Condition on Cost function  if iterative.*/
    PARAM m_minRelSol;/**<Condition on successive estimates if iterative.*/
    PARAM m_gammaProx;/**<Multiplicative factor for function.*/
    PARAM m_precision;/**<Precision.*/
    StopCriteria m_stopping;/**<Choice of stopping condition.*/
    Manifold m_manifold= Manifold::Rn;/**< Manifold flag, only Rn yet.*/
public:
    //constructors/destructors
    /**
      * Standart constructor
    */
    ProxOperator(){init();}

    /**
      * Standart destructor.
    */
    virtual ~ProxOperator(){}

    //Setters
    /**
      * Default initialization: stopping after 10 iterations.
    */
    void init();
    /**
      * Default precision.
    */
    void setDefaultPrecision();

    /**
      * Set Number of iterations
       * @param[in] niter: number of iterations.
    */
    void setNiter(const size_t &niter){m_nIter=niter;}

    /**
      * Set minimal linf norm on gradient
       * @param[in] GradLinf: minimal linf norm on gradient.
    */
    void setMinGradLInf(const PARAM &GradLinf){m_minGradlinf=GradLinf;}

    /**
      * Set minimal relative difference of successive cost function
       * @param[in] MinRelCost: minimal relative cost.
    */
    void setMinRelCost(const PARAM &MinRelCost){m_minRelCost=MinRelCost;}

    /**
      * Set minimal relative difference of successive estimates
       * @param[in] MinRelSol: minimal relative difference in estimates.
    */
    void setMinRelSol(const PARAM &MinRelSol){m_minRelSol=MinRelSol;}

    /**
      * Set multiplication factor for prox (prox_{gamma*f()})
       * @param[in] gamma: multiplication factor for prox .
    */
    void setMultProx(const PARAM& gamma) const {
        if(gamma>0) m_gammaProx=gamma;
        else std::cout<<"Gamma should be >0. Nothing done"<<std::endl;
    }


    /**
      * Set stopping criterion
      * @param[in] Stop: stopping criterion.
    */
    void setStopping(const StopCriteria &Stop){m_stopping=Stop;}


    //Getters

    /**
      * Get Number of iterations
       * @return number of iterations.
    */
    size_t getNiter() const {return m_nIter;}

    /**
      * Get minimal linf norm on gradient
       * @return minimal linf norm on gradient.
    */
    PARAM getMinGradLInf() const {return m_minGradlinf;}

    /**
      * Get minimal relative difference of successive cost function
       * @return minimal relative cost.
    */
    PARAM getMinRelCost() const {return m_minRelCost;}

    /**
      * Get minimal relative difference of successive estimates
       * @return minimal relative difference in estimates.
    */
    PARAM  getMinRelSol() const {return m_minRelSol;}

    /**
      * Get multiplication factor for prox (prox_{gamma*f()})
       * @return multiplication factor for prox .
    */
    PARAM  getMultProx() const {return m_gammaProx;}

    /**
      * Get stopping criterion
      * @return stopping criterion.
    */
    StopCriteria getStopping() const {return m_stopping;}

    /**
      * Interface to get the signature of the object.
       * @return signature of the object.
    */
    virtual std::string getSignature() const =0 ;

    //processing functions
    /**
      * Interface to compute the rprox operator.
       * @param[in] RefData point where to compute the rprox
       * @param[in] Nelems number of elements in reference point
       * @return array containing proximity operator.
    */
    PARAM* computeRProx(PARAM* RefData,const size_t &Nelems);

    /**
      * Interface to compute the prox of a shifted function.
       * @param[in] RefData point where to compute the rprox
       * @param[in] Nelems number of elements in reference point
       * @param[in] CtrRef shift applied to the function
       * @return array containing proximity operator.
    */
    PARAM* computeProxCtr(const PARAM* RefData, const size_t &Nelems,
                                                        const PARAM& CtrRef);

    /**
      * Interface to compute the prox of a weighted function.
       * @param[in] RefData point where to compute the rprox
       * @param[in] Nelems number of elements in reference point
       * @param[in] multfactor factor applied to the function
       * @return array containing proximity operator.
    */
    PARAM* computeProxMult(const PARAM* RefData, const size_t &Nelems,
                                                    const PARAM& multfactor);

    /**
      * Compute Prox Operator of dual function using Moreau identity.
       * @param[in] RefData point where to compute the prox
       * @param[in] Nelems number of elements in reference point
       * @return array containing prox of dual function.
    */
    PARAM* computeProxDualMoreau(const PARAM* RefData,
                                                    const size_t &Nelems);

    /**
      * Compute Prox Operator of dual function.
       * @param[in] RefData point where to compute the prox of  dual
       * @param[in] Nelems number of elements in reference point
       * @return  array containing prox of dual function.
       @note The default version is using the Moreau identity. If known, this
       can be computed directly.
    */
    virtual PARAM* computeProxDual(PARAM* RefData, const size_t &Nelems);

    /**
      * Compute value of function at given point.
       * @param[in] RefData point where to compute the function
       * @param[in] Nelems number of elements in reference point
       * @return  value of function at the given point.
      @note The default version is using the Moreau identity. If known, this
                                                           can be computed directly.
                                                        */
    virtual PARAM computeFunction(PARAM* RefData,const size_t &Nelems)=0;

};

#endif //PROXOP_H

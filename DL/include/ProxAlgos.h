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
 **    Date:  10/2017
 **
 **    File:  ProxAlgos.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces Proximal Algorithms
 **    -----------
 ** - 10/2017: Virtual Class with gradient properties
 ******************************************************************************/

#ifndef PROXALGO_H
#define PROXALGO_H

#include "MatrixDecomposition.h"
#include "S1Decomposition.h"
#include "OpenMPInterface.h"
#include "ProxOperators.h"
#include "GradOperators.h"

#include <limits>
#include <vector>
#include <list>
#include <algorithm>
#include <functional>


template <typename PARAM>  class ProxAlgo {
protected:
    size_t m_nIter;/**<Max number of Iterations*/
    PARAM m_minGradlinf;/**<Gradient condition if iterative.*/
    PARAM m_minRelCost;/**<Condition on Cost function  if iterative.*/
    PARAM m_minRelSol;/**<Condition on successive estimates if iterative.*/
    PARAM m_precision;/**<Precision.*/
    StopCriteria m_stopping;/**<Choice of stopping condition.*/
    bool m_verbose;/**< Verbosity*/
    Manifold m_manifold= Manifold::Rn;/**< Manifold flag, only Rn yet.*/
public:
    /**
      * Standart constructor
    */
    ProxAlgo(GradOperator<PARAM>& GradOp,std::vector<ProxOperator<PARAM> >&
        listProx){
            init();
    }

    /**
      * Standart destructor.
    */
    virtual ~ProxAlgo(){}

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
      * Set stopping criterion
      * @param[in] Stop: stopping criterion.
    */
    void setStopping(const StopCriteria &Stop){m_stopping=Stop;}

    /**
      * Set stopping criterion
      * @param[in] Stop: stopping criterion.
    */
    void setVerbose(const bool _verb){m_verbose=_verb;}


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
      * Get stopping criterion
      * @return stopping criterion.
    */
    StopCriteria getStopping() const {return m_stopping;}

    /**
      * Get stopping criterion
      * @return stopping criterion.
    */
    bool getVerbose() const {return m_verbose;}


    /**
      * Interface to get the signature of the object.
       * @return signature of the object.
    */
    virtual std::string getSignature() const =0 ;


    //processing functions
    /**
      * Interface to minimize the cost functions.
       * @param[in] RefData initialization Point
       * @param[in] Nelems number of elements in initialization point
       * @return array containing solution of minimization problem
    */
    virtual PARAM* minimizeCost(const PARAM* RefData,const size_t &Nelems)=0;


};



/**
 * @class gForwardBackward
 * @brief Generalized Forward-Backward Implementation.
 * @see Raguet et al SIAM 2013
*/
template <typename PARAM>  class gForwardBackward : public ProxAlgo<PARAM>{
protected:
    GradOperator<PARAM> m_gradOp;
    std::vector<ProxOperator<PARAM> >* m_listProx;/**< List of proximity operators.*/
    std::vector<PARAM>* m_weightProx;/**< List of convex weights for each prox.*/
    std::vector<PARAM>* m_lambda;/**< List of relaxation weights across iterations.*/
    std::vector<PARAM>* m_gamma;/**< List of gradient descent steps across iterations.*/
    size_t m_nProx;/**< Number of proximity operators.*/
    bool m_weightOwnFlag;/**< Flag to specify if weights are owned.*/
public:
    //constructors/destructors
    /**
      * Standart constructor
    */
    gForwardBackward(GradOperator<PARAM>* GradOp,std::vector<ProxOperator<PARAM> >*
        listProx){
            init();
            m_weightOwnFlag=false;
            m_listProx=listProx;
            m_gradOp=GradOp;
            m_nProx=listProx.size();
            m_weightProx->assign(m_nProx,1.0/m_nProx);
    }

    /**
      * Standart destructor.
    */
    virtual ~gForwardBackward(){}

    //Setters
    /**
      * Default initialization.
    */
    void init();

    //Getters
    /**
      * Interface to get the signature of the object.
       * @return signature of the object.
    */
    virtual std::string getSignature() {
        std::stringstream strstr;
        std::string Sign=std::string("generalized-Forward-Backward-");
        strstr<<this->m_radius;
        Sign.append(strstr.str());
        return Sign;

    }

    /**
      * Get relaxation parameter for current iteration
      * @param[in] iteration iteration umber
      * @return relaxation parameter.
    */
    inline PARAM& get_lambda(const size_t iteration){
        return m_lambda[MIN(iteration,m_lambda.size())];
    }

    /**
      * Get gradient descent step for current iteration
      * @param[in] iteration iteration umber
      * @return gradient step.
    */
    inline PARAM& get_gamma(const size_t iteration){
        return m_gamma[MIN(iteration,m_gamma.size())];
    }



    //setters
    /**
      * Default precision.
    */
    void setDefaultPrecision();

    /**
      * Set Weights for proximal temporary variables.
    */
    void setWeightProx(std::vector<PARAM>* _weightProx){
        typename std::vector<PARAM>::iterator itPos;
        itPos=std::find_if(_weightProx->begin(),_weightProx->end(),
                                         std::bind2nd(less_equal<PARAM>(),0.0));
        if(itPos!=_weightProx->end()){
            std::string exMessage=EXCEPT_WHERE();
            exMessage.append("Incorrect weight: non positive value: ");
            std::stringstream strstr;
            strstr<<*itPos;
            exMessage.append(strstr.str());
            throw std::invalid_argument(exMessage);
        }

        PARAM total=std::accumulate(_weightProx->begin(),_weightProx->end(),0.0);
        if(_weightProx.size()==m_nProx) {
            if(m_weightOwnFlag) m_weightProx.clear();
            m_weightProx=_weightProx;
            if(fabs(total-1.0)>this->m_precision){
                std::cout<<"Weights for proximity operators not convex. ";
                std::cout<<"Sum:"<<total<<". Divide by sum"<<std::endl;
            }
            itPos=m_weightProx->begin();
            for(size_t kp=0;kp<m_nProx;++kp,++itPos) *itPos/=total;
            m_weightOwnFlag=false;
        } else {
            std::string exMessage=EXCEPT_WHERE();
            exMessage.append("Incorrect weight size: ");
            std::stringstream strstr;
            strstr<<_weightProx.size()<<" vs "<<m_nProx;
            exMessage.append(strstr.str());
            throw std::invalid_argument(exMessage);
        }
    }

    /**
      * Set iteration dependent weights for proximal temporary variables.
    */
    void setRelaxWeights(std::vector<PARAM>* _lambda){ m_lambda=_lambda;}

    /**
      * Set iteration dependent grad steps for proximal temporary variables.
    */
    void setGradSteps(std::vector<PARAM>* _gamma){ m_gamma=_gamma;}

    //processing functions
    /**
      * Interface to minimize the cost functions.
       * @param[in] RefData initialization Point
       * @param[in] Nelems number of elements in initialization point
       * @return array containing solution of minimization problem
    */
    virtual PARAM* minimizeCost(const PARAM* RefData,const size_t &Nelems){
        //To code
    };

};


#endif //PROXALGO_H

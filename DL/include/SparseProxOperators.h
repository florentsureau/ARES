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
 **    File:  SparseProxOperators.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces Proximity Operator for sparsity enforcing functions
 **    -----------
 ** - 09/2017: L1, "L0" with weighted counterparts, L1+affine, L1+simplex
 ******************************************************************************/

 #ifndef SPARSEPROXOP_H
 #define SPARSEPROXOP_H

 #include "MatrixDecomposition.h"
 #include "S1Decomposition.h"
 #include "OpenMPInterface.h"

 #include <limits>
 #include <vector>
 #include <list>
 #include <algorithm>
 #include <functional>


 #include "ProxOperators.h"

/**
 * @class ProxL1
 * @brief Class implementing proximity operator for L1 norm. The weighting is
 * given by @see ProxOperator::m_GammaProx

*/
template <typename PARAM>  class ProxL1:public  ProxOperator<PARAM> {
protected:
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    ProxL1(){}
    /**
      * Standart destructor.
    */
    virtual ~ProxL1();

    /**
      * Signature of ProxL1
       * @return signature of the object.
    */
    virtual std::string getSignature(){
        std::stringstream strstr;
        std::string Sign=std::string("Prox-L1-gamma");
        strstr<<ProxOperator<PARAM>::getMultProx();
        Sign.append(strstr.str());
        return Sign;
    };

    /**
      * Interface to compute the proximity operator of l1 norm - soft thresholding
        * @param[in] RefData point where to compute the prox
        * @param[in] Nelems number of elements in reference point
        * @return array containi.
    */
    virtual PARAM* computeProx(PARAM* data,const size_t &Nelems);

    /**
      * Compute Prox Operator of dual: projection onto linf ball
       * @return  array containing prox of dual function.
       @note This is the non-positive projection.
    */
    virtual PARAM*  computeProxDual(PARAM* RefData, const size_t &Nelems);

    /**
      * Compute value of function at given point.
       * @param[in] RefData point where to compute the function
       * @param[in] Nelems number of elements in reference point
       * @return  value of function at the given point.
      @note The default version is using the Moreau identity. If known, this
                                                           can be computed directly.
                                                        */
    virtual PARAM computeFunction(PARAM* RefData,const size_t &Nelems);
};

/**
 * @class ProxL2
 * @brief Class implementing proximity operator for L2 norm. The weighting is
 * given by @see ProxOperator::m_GammaProx

*/
template <typename PARAM>  class ProxL2:public  ProxOperator<PARAM> {
protected:
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    ProxL2(){}
    /**
      * Standart destructor.
    */
    virtual ~ProxL2();

    /**
      * Signature of ProxL2
       * @return signature of the object.
    */
    virtual std::string getSignature(){
        std::stringstream strstr;
        std::string Sign=std::string("Prox-L2-gamma");
        strstr<<ProxOperator<PARAM>::getMultProx();
        Sign.append(strstr.str());
        return Sign;
    };

    /**
      * Interface to compute the proximity operator of l2 norm
        * @param[in] RefData point where to compute the prox
        * @param[in] Nelems number of elements in reference point
        * @return array containi.
    */
    virtual PARAM* computeProx(PARAM* data,const size_t &Nelems);

    /**
      * Compute Prox Operator of dual: projection onto l2 ball
       * @return  array containing prox of dual function.
       @note This is the non-positive projection.
    */
    virtual PARAM*  computeProxDual(PARAM* RefData, const size_t &Nelems);

    /**
      * Compute value of function at given point.
       * @param[in] RefData point where to compute the function
       * @param[in] Nelems number of elements in reference point
       * @return  value of function at the given point.
      @note The default version is using the Moreau identity. If known, this
                                                           can be computed directly.
                                                        */
    virtual PARAM computeFunction(PARAM* RefData,const size_t &Nelems);
};


/**
 * @class ProjSparseSimplex
 * @brief Class implementing sparse projection onto the simplex.
*/
template <typename PARAM>  class ProjSparseSimplex:public  ProxOperator<PARAM> {
protected:
    size_t m_sparsity;
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    ProjSparseSimplex():m_sparsity(1){}
    /**
      * Standart destructor.
    */
    virtual ~ProjSparseSimplex();

    /**
      * Set sparsity level
       * @param[in] ksparse sparsity level.
    */
    void setSparsity(const PARAM & ksparse){
        if(ksparse>0) m_sparsity=ksparse;
        else std::cout<<"Sparsity should be >0. Nothing done"<<std::endl;
    }

    /**
      * Signature of ProjSparseSimplex
       * @return signature of the object.
    */
    virtual std::string getSignature(){
        std::stringstream strstr;
        std::string Sign=std::string("Proj-Sparse-Simplex-k");
        strstr<<m_sparsity;
        Sign.append(strstr.str());
        return Sign;};


    /**
      * Get sparsity level
       * @return sparsity level.
    */
    void getSparsity() const {return m_sparsity;}


    /**
      * Interface to compute the sparse projection onto the simplex.
       * @return array containing sparse projection onto the simplex.
    */
    virtual PARAM* computeProx(PARAM* data,const size_t &Nelems);

    /**
      * Compute Prox Operator of dual of the sparse projection onto the simplex.
       * @return  array containing prox of dual function.
       @note This is using the Moreau Identity.
    */
    virtual PARAM* computeProxDual(PARAM* RefData, const size_t &Nelems){
        return ProxOperator<PARAM>::computeProxDual(RefData,Nelems);
    }

    /**
      * Compute value of function at given point.
       * @param[in] RefData point where to compute the function
       * @param[in] Nelems number of elements in reference point
       * @return  value of function at the given point.
      @note The default version is using the Moreau identity. If known, this
                                                           can be computed directly.
                                                        */
    virtual PARAM computeFunction(PARAM* RefData,const size_t &Nelems);
};

/**
 * @class ProjSparseAffine
 * @brief Class implementing sparse projection on affine constraints.
*/
template <typename PARAM>  class ProjSparseAffine:public  ProxOperator<PARAM> {
protected:
    size_t m_sparsity;
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    ProjSparseAffine():m_sparsity(1){}
    /**
      * Standart destructor.
    */
    virtual ~ProjSparseAffine();

    /**
      * Set sparsity level
       * @param[in] ksparse sparsity level.
    */
    void setSparsity(const PARAM & ksparse){
        if(ksparse>0) m_sparsity=ksparse;
        else std::cout<<"Sparsity should be >0. Nothing done"<<std::endl;
    }

    /**
      * Signature of ProjSparseAffine
       * @return signature of the object.
    */
    virtual std::string getSignature(){
        std::stringstream strstr;
        std::string Sign=std::string("Proj-Sparse-Simplex-k");
        strstr<<m_sparsity;
        Sign.append(strstr.str());
        return Sign;};


    /**
      * Get sparsity level
       * @return sparsity level.
    */
    void getSparsity() const {return m_sparsity;}


    /**
      * Interface to compute the sparse projection onto the simplex.
       * @return array containing sparse projection onto the simplex.
    */
    virtual PARAM* computeProx(PARAM* data,const size_t &Nelems);

    /**
      * Compute Prox Operator of dual of the sparse projection onto the simplex.
       * @return  array containing prox of dual function.
       @note This is using the Moreau Identity.
    */
    virtual PARAM* computeProxDual(PARAM* RefData, const size_t &Nelems){
        return ProxOperator<PARAM>::computeProxDual(RefData,Nelems);
    }

    /**
      * Compute value of function at given point.
       * @param[in] RefData point where to compute the function
       * @param[in] Nelems number of elements in reference point
       * @return  value of function at the given point.
      @note The default version is using the Moreau identity. If known, this
                                                           can be computed directly.
                                                        */
    virtual PARAM computeFunction(PARAM* RefData,const size_t &Nelems);
};

#endif //SPARSEPROXOP_H

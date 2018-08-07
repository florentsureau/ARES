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
 **    File:  ProjOperators.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces Proximity Operator of indicator functions of
 **    convex sets (projections)
 **    -----------
 ** - 09/2017: Virtual Class
 ******************************************************************************/

 #include "ProxOperators.h"

/**
 * @class ProjNonNegOrthant
 * @brief Class implementing the non-negative projection.
*/
template <typename PARAM>  class ProjNonNegOrthant:public  ProxOperator<PARAM> {
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    ProjNonNegOrthant(){}
    /**
      * Standart destructor.
    */
    virtual ~ProjNonNegOrthant();

    /**
      * Signature of ProjNonNegOrthant
       * @return signature of the object.
    */
    virtual std::string getSignature(){
        std::string Sign=std::string("Proj-NonNegative-Orthant");
        return Sign;
    };

    /**
      * Interface to compute the proximity operator of the indicator of the
        non-negative orthant.
        * @param[in] RefData point where to compute the prox
        * @param[in] Nelems number of elements in reference point
        * @return array containi.
    */
    virtual PARAM* computeProx(PARAM* data, const size_t &Nelems);

    /**
      * Compute Prox Operator of dual non-negative projection.
       * @param[in] RefData point where to compute the prox of dual
       * @param[in] Nelems number of elements in reference point
       * @return  array containing prox of dual function.
       @note This is the non-positive projection.
    */
    virtual PARAM* computeProxDual(PARAM* RefData,const size_t &Nelems);

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
 * @class ProjNonPosOrthant
 * @brief Class implementing non-negative projection.
*/
template <typename PARAM>  class ProjNonPosOrthant:public  ProxOperator<PARAM>{
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    ProjNonPosOrthant(){}
    /**
      * Standart destructor.
    */
    virtual ~ProjNonPosOrthant();

    /**
      * Signature of ProjNonPosOrthant
       * @return signature of the object.
    */
    virtual std::string getSignature(){
        std::string Sign=std::string("Proj-NonPositive-Orthant");
        return Sign;
    };

    /**
      * Interface to compute the proximity operator of the indicator of
        the non-positive orthant.
        * @param[in] RefData point where to compute the prox
        * @param[in] Nelems number of elements in reference point
       * @return array containing projection onto the non-positive orthant.
    */
    virtual PARAM* computeProx(PARAM* data,const size_t &Nelems);

    /**
      * Compute Prox Operator of dual of non-positive projection.
       * @param[in] RefData point where to compute the prox of dual
       * @param[in] Nelems number of elements in reference point
       * @return  array containing prox of dual function.
       @note This is the non-positive projection.
    */
    virtual PARAM* computeProxDual(PARAM* RefData,const size_t &Nelems);

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
 * @class ProjL2Ball
 * @brief Class implementing projection onto the L2 Ball with radius(gamma).
*/
template <typename PARAM>  class ProjL2Ball:public ProxOperator<PARAM> {
protected:
    PARAM m_radius;
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    ProjL2Ball():m_radius(1){}
    /**
      * Standart destructor.
    */
    virtual ~ProjL2Ball();

    /**
      * Set radius of L2 ball
       * @param[in] Radius: radius of the ball.
    */
    void setRadius(const PARAM & Radius){
        if(m_radius>0) m_radius=Radius;
        else std::cout<<"Radius should be >0. Nothing done"<<std::endl;
    }

    /**
      * Signature of ProjL2Ball
       * @return signature of the object.
    */
    virtual std::string getSignature(){
        std::stringstream strstr;
        std::string Sign=std::string("Proj-L2Ball-Radius");
        strstr<<m_radius;
        Sign.append(strstr.str());
        return Sign;
    };

    /**
      * Get radius of L2 ball
       * @return radius of the ball.
    */
    void getRadius() const {return m_radius;}

    /**
      * Interface to compute the proximity operator of the l2 ball constraint.
       * @param[in] RefData point where to compute the prox
       * @param[in] Nelems number of elements in reference point
       * @return array containing projection onto the l1 ball.
    */
    virtual PARAM* computeProx(PARAM* data, const size_t &Nelems);

    /**
      * Compute Prox Operator of dual of l2 ball.
       * @param[in] RefData point where to compute the prox of dual
       * @param[in] Nelems number of elements in reference point
       * @return  array containing prox of dual function.
      @note This is using the Moreau Identity.
    */
    virtual PARAM* computeProxDual(PARAM* RefData,const  size_t &Nelems){
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
 * @class ProjLInfBall
 * @brief Class implementing projection onto the LInf Ball with radius(gamma).
*/
template <typename PARAM>  class ProjLInfBall:public ProxOperator<PARAM> {
protected:
    PARAM m_radius;
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    ProjLInfBall():m_radius(1){}
    /**
      * Standart destructor.
    */
    virtual ~ProjLInfBall();

    /**
      * Set radius of L2 ball
       * @param[in] Radius: radius of the ball.
    */
    void setRadius(const PARAM & Radius){
        if(m_radius>0) m_radius=Radius;
        else std::cout<<"Radius should be >0. Nothing done"<<std::endl;
    }

    /**
      * Signature of ProjLInfBall
       * @return signature of the object.
    */
    virtual std::string getSignature(){
        std::stringstream strstr;
        std::string Sign=std::string("Proj-LInfBall-Radius");
        strstr<<m_radius;
        Sign.append(strstr.str());
        return Sign;
    };

    /**
      * Get radius of Linf ball
       * @return radius of the ball.
    */
    void getRadius() const {return m_radius;}

    /**
      * Interface to compute the proximity operator of the linf ball constraint.
       * @param[in] RefData point where to compute the prox
       * @param[in] Nelems number of elements in reference point
       * @return array containing projection onto the l1 ball.
    */
    virtual PARAM* computeProx(PARAM* data,const size_t &Nelems);

    /**
      * Compute Prox Operator of dual of linf ball.
       * @param[in] RefData point where to compute the prox of dual
       * @param[in] Nelems number of elements in reference point
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
 * @class ProjAffine
 * @brief Class implementing projection onto the affine constraint.
*/
template <typename PARAM>  class ProjAffine:public  ProxOperator<PARAM> {
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    ProjAffine(){}
    /**
      * Standart destructor.
    */
    virtual ~ProjAffine();

    /**
      * Signature of ProjAffine
       * @return signature of the object.
    */
    virtual std::string getSignature(){
        std::string Sign=std::string("Proj-Affine");
        return Sign;
    };

    /**
      * Interface to compute the proximity operator of the affine constraint.
       * @param[in] RefData point where to compute the prox
       * @param[in] Nelems number of elements in reference point
       * @return array containing projection onto the affine constraint.
    */
    virtual PARAM* computeProx(PARAM* data, const size_t &Nelems);

    /**
      * Compute Prox Operator of dual of affine constraint.
       * @param[in] RefData point where to compute the prox of dual
       * @param[in] Nelems number of elements in reference point
       * @return  array containing prox of dual function.
      @note This is using the Moreau Identity.
    */
    virtual PARAM* computeProxDual(PARAM* RefData,const size_t &Nelems){
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
 * @class ProjSimplex
 * @brief Class implementing projection onto the simplex.
*/
template <typename PARAM>  class ProjSimplex: public ProxOperator<PARAM> {
public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    ProjSimplex(){}
    /**
      * Standart destructor.
    */
    virtual ~ProjSimplex();

    /**
      * Signature of ProjSimplex
       * @return signature of the object.
    */
    virtual std::string getSignature(){
        std::string Sign=std::string("Proj-Simplex");
        return Sign;
    };

    /**
      * Interface to compute the proximity operator of the indicator of
        the simplex.
        * @param[in] RefData point where to compute the prox
        * @param[in] Nelems number of elements in reference point
        * @return array containing projection onto the simplex.
    */
    virtual PARAM* computeProx(PARAM* data, const size_t &Nelems);

    /**
      * Compute Prox Operator of dual of  indicator of
        the simplex.
        * @param[in] RefData point where to compute the prox of dual
        * @param[in] Nelems number of elements in reference point
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

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
 **    File:  DictionaryInit.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces classes for Dictionary Initialization.
 **    -----------
 **
 **
 ******************************************************************************/


#ifndef DictionaryInit_H
#define DictionaryInit_H
#include "TempArray.h"
#include "MatrixOper.h"
#include "IM_IO.h"
#include "gslRNG.h"

enum class InitDictionaryType {
    RandomAtomsDictionary,
    RandomSparseDictionary,
    UserDefinedDictionary,
};

class DictionaryInit {
public:
    //Constructors/Destructors
    DictionaryInit(const size_t Nx,const size_t Ny){
        m_Nx = Nx;
        m_Ny = Ny;
    }
    virtual ~DictionaryInit(){}
    //Setters
    void setRNG(gslRNG* RNG) {p_RNG = RNG;}
    //Processing functions
    virtual dblarray* computeInit()=0;
    virtual std::string getSignature() const =0;
    //IOs
    static dblarray* readDicoFits(char *FitsFileName){
        dblarray* InitDict=new dblarray();
        fits_read_dblarr (FitsFileName, *InitDict);
        return InitDict;}
    static void writeDicoFits(char *FitsFileName,dblarray Dico){
        fits_write_dblarr(FitsFileName, Dico);}
protected:
    size_t m_Nx;
    size_t m_Ny;
    gslRNG* p_RNG;
private:
};


class SparseDictInit: public DictionaryInit {
public:
    //Constructors/Destructors
    SparseDictInit(const double density,const size_t Nx,const size_t Ny):
        DictionaryInit(Nx,Ny){
        m_density=density;
        m_Nx=Nx;
        m_Ny=Ny;
    }
    virtual ~SparseDictInit(){}
    void initRNG(){
        if(rng !=NULL) gsl_rng_free (rng);
        const gsl_rng_type * T;
        T = gsl_rng_default;
        rng = gsl_rng_alloc (T);
    }
    virtual std::string getSignature() const {
        return std::string("SparseInit");
    }

    //Processing functions
    virtual dblarray* computeInit();
protected:
    gsl_rng* rng;
    double m_density;
};

class RandomDictInit: public DictionaryInit {
public:
    //Constructors/Destructors
    RandomDictInit(const size_t Natoms,dblarray* Training);
    virtual ~RandomDictInit() {delete m_training;}
    void initRNG() {
        if(rng !=NULL) gsl_rng_free (rng);
        const gsl_rng_type * T;
        T = gsl_rng_default;
        rng = gsl_rng_alloc (T);
    }
    virtual std::string getSignature() const {
        return std::string("RandomTrainInit");
    }
    //Processing functions
    virtual dblarray* computeInit();
protected:
    gsl_rng* rng;
    double m_density;
    dblarray* m_training;
};
#endif
//DictionaryInit_H

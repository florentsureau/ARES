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
 **    File:  DictionaryInit.cpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Implementation of dictionary learning initialization.
 **    -----------
 **  - 07/2016: Sparse Unif. Rand. Initialization, reading routines
 **
 ******************************************************************************/

#include "DictionaryInit.h"
#include "myexit.hpp"

RandomDictInit::RandomDictInit(const size_t Natoms, dblarray* Training) : DictionaryInit(0,0) {
    std::cout<<"TRAINING"<<std::endl;
    m_Nx = Natoms;
    m_Ny = Training->nx();
    m_training = new dblarray();
    m_training->alloc(Training->buffer(),Training->nx(),Training->ny(),Training->nz());
}

dblarray* SparseDictInit::computeInit(){
    size_t Nelem= m_Nx* m_Ny,Nkeep;
    if((m_density >=1)||(m_density<=0)) {
        std::cout<<"Sparse density should be 0<d<1:"<<m_density <<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    gsl_permutation* Perm=p_RNG->createPermutation(Nelem);
    Nkeep=(int) Nelem* m_density;

    dblarray* InitRandom=new dblarray(m_Nx, m_Ny);
    double *buffer=InitRandom->buffer();
    InitRandom->init(0.);
    for(size_t k=0;k<Nkeep;++k)
        buffer[gsl_permutation_get(Perm,k)]= (double) p_RNG->getBernoulli(0.5)-0.5;
    return InitRandom;
}

dblarray* RandomDictInit::computeInit(){
    if((size_t)m_training->ny() < m_Nx) {
        std::cerr << "/!\\/!\\ Error : the size of the training set (" << m_training->ny() << ") must be larger than the number of dictionary atoms (" << m_Nx << ")" << std::endl;
        return new dblarray(0);
    }
    gsl_permutation* Perm=p_RNG->createPermutation(m_training->ny());
    dblarray* InitRandom=new dblarray(m_Nx, m_Ny);
    double *bufferInit, *bufferTrain;
    for(size_t k=0;k<m_Nx;++k) {
        bufferTrain = m_training->buffer()+gsl_permutation_get(Perm,k)*m_Ny;
        bufferInit= InitRandom->buffer()+k;
        for(size_t ky=0;ky<m_Ny;++ky, bufferInit +=m_Nx,++bufferTrain)
            *bufferInit= *bufferTrain;
    }
    for(size_t k=0;k<m_Nx;++k)
        std::cout<<"Atom "<<k<<"="<<"Train "<<gsl_permutation_get(Perm,k)<<std::endl;
    return InitRandom;
}

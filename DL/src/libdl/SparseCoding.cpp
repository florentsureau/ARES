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
 **    File:  SparseCoding.cpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Implementation of abstract class for sparse decomposition
 **    -----------
 **
 ******************************************************************************/


#include "SparseCoding.h"
#include <vector>
#include <omp.h>

#include "IM_IO.h"
#include "IM1D_IO.h"
#include "ISAP_verbose.hpp"
#include "myexit.hpp"

#include "TimeInfoManager.h"

#if ENABLE_VERBOSE_DEBUG_IN_SparseCoding
    #define DEBUGME_(x) x
#else
    #define DEBUGME_(x)
#endif

// ============================
//     SparseCoding_Settings
// ============================

SparseCoding_Settings::SparseCoding_Settings() {
    init();
}
void SparseCoding_Settings::init() {
    verbose           = false;
    errorTarget       = 0;
    sparsityTarget    = 0;
    minimalSparsity   = 0;
    minimalCorrelation= 0.f;
    meanSparsityTarget= 0.f;
    epsilon           = EpsDoublePrec;
    m_manifold        = Manifold::Rn;
}
void SparseCoding_Settings::print(std::string sparseCodingName) const {
    std::cout << "SparseCoding_Settings (" << sparseCodingName << ")" << std::endl;
    std::cout << " errorTarget          = " << errorTarget     << std::endl;
    std::cout << " sparsityTarget       = " << sparsityTarget  << std::endl;
    std::cout << " minimalSparsity      = " << minimalSparsity << std::endl;
    std::cout << " minimalCorrelation   = " << minimalCorrelation << std::endl;
    std::cout << " meanSparsityTarget   = " << meanSparsityTarget << std::endl;
    std::cout << " epsilon              = " << epsilon         << std::endl;
}


// ===================
//     SparseCoding
// ===================
SparseCoding::SparseCoding() {
    init();
}
SparseCoding::SparseCoding(dblarray &dictionary){
    init();
    setDictionary(dictionary);
}
SparseCoding::~SparseCoding(){
    if(m_dictionary!= nullptr) {
        DEBUG_( printf("~SparseCoding() : delete dictionary\n"); );
        delete m_dictionary;
    }
    m_dictionary= nullptr;
    if(true == m_matdecompIsLocal)
        delete m_matdecomp;
}
void SparseCoding::init() {
    m_dictionary= nullptr;
    m_matdecomp=nullptr;
    m_natoms=0;
    m_npix=0;
    m_ninput=0;
    m_boolInputMetric=false;
    m_boolDiagInputMetric=false;
    m_boolCodeMetric=false;
    m_boolDiagCodeMetric=false;
    m_settings.init();
}
void SparseCoding::print() {
    m_settings.print(getName());
}
//i/o
dblarray* SparseCoding::readMetric(const char *FitsFileName){
    dblarray* InitMetric=new dblarray();
    fits_read_dblarr (FitsFileName, *InitMetric);
    return InitMetric;
}

//setters
void SparseCoding::setDictionary(const dblarray &dictionary) {
    if(m_dictionary!=nullptr) delete m_dictionary;
    m_dictionary=new dblarray(dictionary);
    size_t nx = m_dictionary->nx();
    if(m_settings.sparsityTarget <= 0 || m_settings.sparsityTarget > nx)
        m_settings.sparsityTarget = nx;
    m_natoms = dictionary.nx();
    m_npix = dictionary.ny();
}
void SparseCoding::setRFactor(const dblarray &rescale_factor) {
    m_rfactor = rescale_factor;
}

int SparseCoding::checkAllocation() {
    if(~m_isAllocated) {
        if(m_dictionary==nullptr || m_matdecomp==nullptr) {
            ERROR_(std::cerr << "Error: SparseCoding is not properly allocated."
                   <<" Check dictionary/matrixDecomposition" << std::endl; );
            return EXIT_FAILURE;
        }
        else
            m_isAllocated = true;
    }
    return EXIT_SUCCESS;
}

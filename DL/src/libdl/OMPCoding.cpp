/******************************************************************************
 **                   Copyright (C) 2017 by CEA
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
 **    File:  OMPCoding.cpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Implementation of Sparse decomposition
 **    -----------
 ** - 06/2016: Implementation of OMP
 ** - 01-02/2017: OMP for S1 data added
 ** - 06/2017: Product Space Sparse coding
 ******************************************************************************/

#include "OMPCoding.h"
//#define DEBUGME_(x) x
#define DEBUGME_(x)
#include <unistd.h>


// ===================
//     OMP
// ===================

OMP::~OMP() {
    if((m_weightedDictionary!=nullptr)&&(m_weightedDictionary!=m_dictionary))
                                                    delete m_weightedDictionary;
    m_weightedDictionary=nullptr;
}

void OMP::setDictionary(const dblarray &dictionary) {
    FINFO_;
    DEBUGME_(
    std::cout<<"SET DICO IN OMP"<<std::endl;);
    SparseCoding::setDictionary(dictionary);
    m_weightedDictionary=m_dictionary;
    m_natoms = m_dictionary->nx();
    m_npix = m_dictionary->ny();
    m_rfactor.alloc(m_natoms);
    m_rfactor.init(1.0);
    DEBUGME_(
    std::cout << " OMP:setDictionary : natoms = " <<
            m_natoms << ", atomLength = " << m_npix << std::endl;);
    FOFNI_;
}

void OMP::setRFactor(const dblarray &rescale_factor) {

    if((size_t) rescale_factor.nx() != m_natoms){
        std::cout << "OMP::Wrong size of rescale_factor "<< rescale_factor.nx()
                  << " vs " << m_natoms <<std::endl;
        std::cout << "Do not Use rescaling" <<std::endl;

    } else SparseCoding::setRFactor(rescale_factor);
}

void OMP::setInputMetric(const dblarray &metric) {
    if(metric.naxis()==1) { //Diagonal case
        std::cout<<"Diagonal Metric for OMP"<<std::endl;
        if ((size_t)metric.nx() != m_npix) {
            std::cout << "OMP::Wrong size of Diagonal Metric: "<< metric.nx()
                      << " vs " << m_npix <<std::endl;
            MyExit::exit(EXIT_FAILURE);
        } else {
            m_inputMetric=metric;
            if(m_inputMetric.min()<=0) {
                std::cout << "OMP:Diag Metric should be strictly positive: "<<
                             m_inputMetric.min()<<std::endl;
                MyExit::exit(EXIT_FAILURE);
            }
            m_boolInputMetric=true;
            m_boolDiagInputMetric=true;
        }
    } else if(metric.naxis()==2){
        std::cout<<"2D Metric for OMP"<<std::endl;
        m_boolInputMetric=true;
        if(((size_t)metric.nx()!=m_npix)||((size_t)metric.ny()!=m_npix)) {
            std::cout << "OMP::Wrong size of Metric: "<< metric.nx()<<"x"<<
                                    metric.ny()<< " vs " << m_npix <<std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        m_inputMetric=metric;
        m_boolDiagInputMetric=false;
    } else {
        std::cout << "OMP::Metric should be of dim 1 or 2 "<< metric.nx()
                  << " vs " << m_npix <<". Metric not used"<<std::endl;
        m_boolInputMetric=false;
    }
}

void OMP::setCodeMetric(const dblarray &metric) {
    if(metric.naxis()==1) { //Diagonal case
        if ((size_t)metric.nx() != m_npix) {
            std::cout << "OMP::Wrong size of Diagonal Metric: "<< metric.nx()
                      << " vs " << m_npix <<std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        m_boolCodeMetric=true;
        m_boolDiagCodeMetric=true;
        m_codeMetric=metric;
        if(m_codeMetric.min()<=0) {
            std::cout << "OMP:Diag Metric should be strictly positive: "<<
                         m_codeMetric.min()<<std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
    } else if(metric.naxis()==2){
        if (((size_t)metric.nx() != m_npix)||((size_t)metric.ny() != m_npix)) {
            std::cout << "OMP::Wrong size of 2D Metric: "<< metric.nx()
                      << metric.ny() << " vs " << m_npix <<std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        m_boolCodeMetric=true;
        m_boolDiagCodeMetric =false;
        m_codeMetric=metric;//copy constructor
        //Cholesky factorization ? Ensure posdef ?
    } else {
        m_boolCodeMetric=false;
        std::cout << "OMP::Metric should be of dim 1 or 2 "<< metric.nx()
                  << " vs " << m_npix <<". Metric not used"<<std::endl;
    }
}

void OMP::setMatDecomp(MatrixDecomposition &mdecomp) {
    FINFO_;
    if((m_matdecomp!=NULL)&&(m_matdecompIsLocal)) {
        DEBUGME_( printf("DELETE DICTIONARY\n"););
        delete m_matdecomp;
    }
    m_matdecomp =&(mdecomp);
    m_matdecompIsLocal = false;
    FOFNI_;
}

void OMP::setMatDecomp() {
    FINFO_;
#if defined _ARMA_DEC
    m_matdecomp = new ArmadilloDecomposition();
#elif defined _ATLAS_DEC
    m_matdecomp = new AtlasDecomposition();
#else
    m_matdecomp = new SapDecomposition();
#endif
    m_matdecompIsLocal = true;
    FOFNI_;
}

std::string OMP::getSignature() {
    return std::string("OMP with decomposition ") + m_matdecomp->getSignature();
}

double OMP::inputDotProduct(const dblarray &vector1, const dblarray &vector2){
    return m_matdecomp->dotProduct(vector1,vector2,m_inputMetric,
                                   m_boolInputMetric,m_boolDiagInputMetric);
}

double OMP::atomsDotProduct(double *vec1, double *atom, double *inmetric){
    return m_matdecomp->dotProduct(vec1, atom, m_natoms,1, m_npix, inmetric,
                                   m_boolInputMetric,m_boolDiagInputMetric);
}


double OMP::codeDotProduct(const dblarray &vector1, const dblarray &vector2){
    return m_matdecomp->dotProduct(vector1,vector2,m_codeMetric,
                                   m_boolCodeMetric,m_boolDiagCodeMetric);
}

double OMP::inputNorm(const dblarray &vector1){
    return sqrt(inputDotProduct(vector1, vector1));
}

double OMP::codeNorm(const dblarray &vector1){
    return sqrt(codeDotProduct(vector1, vector1));
}

int OMP::findClosestAtom(double* residualPtr,double* codePtr,
                                      double *inputMetricPtr, double* scalePtr){
    // Find maximum correlating atom
    int new_pos; //new index in support
    double* iterPtr = m_dictionary->buffer();
    for(unsigned int k=0; k<m_natoms;++k,++iterPtr) {
        codePtr[k] = atomsDotProduct(iterPtr, residualPtr, inputMetricPtr);
        //codeptr[k]= m_matdecomp->dotProduct(iterptr,resptr,m_natoms,1,
        //          m_npix,inmetricptr,m_boolinputmetric,m_booldiaginputmetric);
    }
    double maxval=0;
    new_pos=0;
    for (size_t k=0; k<m_natoms;++k) {
        codePtr[k] /= scalePtr[k];
        if(fabs(codePtr[k])>maxval){
            new_pos=k;
            maxval=fabs(codePtr[k]);
        }
    }
    return(new_pos);
}

void OMP::getInputSample(DataStruc*& Input,size_t sampleNb,
                                                        double *InSample){
    size_t Npatches = Input->getNPatches();
    size_t dimChannel=Input->getManifoldDimSize(m_channel);
    if((Input==NULL)||(sampleNb> Npatches)){
        std::cout<<"Wrong input sample required"<<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    double *iterBuf=Input->getSampleChannel(sampleNb,m_channel);
    double *iterIn= InSample;
    for(size_t np=0;np<dimChannel;++np,++iterBuf,++iterIn) *iterIn=*iterBuf;
}

dblarray* OMP::applySqrtWeights(dblarray*& Input, dblarray& Weights){
    size_t pixNb=Input->nx();
    if(pixNb!=m_npix){
        std::cout<<"Wrong input sample required"<<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    if(m_dictionary==nullptr) {
        std::cout<<"Should first set the dictionary"<<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    if((m_weightedDictionary!=nullptr)&&(m_weightedDictionary!=m_dictionary))
                                                    delete m_weightedDictionary;
    dblarray* wInput=new dblarray(*Input);
    if(Weights.naxis()==1) { //Diagonal case
        std::cout<<"m_boolDiagInputMetric"<<std::endl;
        m_weightedDictionary=new dblarray(*m_dictionary);
        dblarray TempSqrtDiag=Weights;
        double * TempSqrtDiagBuf=TempSqrtDiag.buffer();
        for(size_t kx=0;kx<(size_t) TempSqrtDiag.n_elem();++kx,++TempSqrtDiagBuf)
            *TempSqrtDiagBuf=sqrt(*TempSqrtDiagBuf);
        std::cout<<"TempSqrtDiag[100]"<<(TempSqrtDiag)(100)<<std::endl;
        m_matdecomp->multRow2DData(*m_weightedDictionary,TempSqrtDiag);
        m_matdecomp->multCol2DData(*wInput,TempSqrtDiag);
    } else if(Weights.naxis()==2) {{
            std::cout<<"m_boolInputMetric"<<std::endl;
            dblarray* TempSqrtDiag=m_matdecomp->sqrtSVD(Weights);
            std::cout<<"TempSqrtDiag[100]"<<(*TempSqrtDiag)(100)<<std::endl;
            m_weightedDictionary=new dblarray(*m_dictionary);
            m_matdecomp->matMult(*m_weightedDictionary,*TempSqrtDiag,*m_dictionary);
            m_matdecomp->matMult(*wInput,*Input,*TempSqrtDiag);
            delete TempSqrtDiag;
        }
    } else {
        std::cout<<"No Preprocessing"<<std::endl;
        m_weightedDictionary=m_dictionary;
    }
    return wInput;
}



/*!
 * \brief OMP::computeCode
 * \param inputVector contains the array (Npix, m_ninput) to sparse code
 * \return the sparse code (m_natoms, m_ninput) of inputVector
 */
dblarray* OMP::computeCode(DataStruc* InputData){
    FINFO_;
    if(checkAllocation()==EXIT_FAILURE) return nullptr;
    if(getManifold()!= InputData->getManifold(m_channel)){
        std::cout<<"BEWARE: INCONSITENCY BETWEEN CODER MANIFOLD "<<getManifold()
                <<" AND INPUT SET MANIFOLD "<< InputData->getManifold(m_channel)
                                                                    <<std::endl;
    }
    //For efficient computation, inputVector should be dblarray(Npix, m_ninput)
    // Code is dblarray(Nindex)
    // Resptr is dblarray(Npix)
    // subDico is dblarray(Nindex,Npix)
    // PinvSubDico is dblarray(Npix, Nindex)
    double *residualPtr=nullptr, *codePtr=nullptr;
    double *inputMetricPtr=nullptr,*iterPtr=nullptr, *vectorPtr=nullptr;
    double *scalePtr = m_rfactor.buffer();
    m_ninput = InputData->getNPatches();
    size_t dimChannel=InputData->getManifoldDimSize(m_channel);
    dblarray* allCodes = new dblarray(m_natoms, m_ninput);//contains all codes
    allCodes->init(0.);
    double residualL2Norm;
    double totalSparsity = 0, processedVec=0;
    int new_pos; //new index in support
    bool stopPursuit=false;

    if (m_npix != (size_t) dimChannel) {
        cout << "Error:Input vector and dictionary dimensions do not match. "
             <<std::endl;
        cout << "   Dico.nx = " << m_npix << ", InputVector->nx = " <<
                dimChannel <<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    if (m_ninput <1){
        std::cout << "OMP::computeCode input vector 2nd dim: " <<
                     m_ninput << " < 1"<< std::endl;
        MyExit::exit(EXIT_FAILURE);
    } else std::cout<<"Process Ninputs="<<m_ninput<<std::endl;
    if (m_settings.verbose == true)
        std::cout << "Running OMP with SparsityTarget = "<<
                        m_settings.sparsityTarget << " and ErrorTarget = " <<
                                            m_settings.errorTarget << std::endl;

    if(m_boolInputMetric) inputMetricPtr = m_inputMetric.buffer();
    else inputMetricPtr =NULL;
    //Weight by Weight^1/2 the dictionary and inputs
    dblarray* inputVectors;
    inputVectors=InputData->getManifoldData(m_channel);

    //Variables used in threading
    size_t Nthreads=OMPThreadHandling::getNThreads();
    char tempname[Nthreads][1024];
    dblarray **pinvsubDico=new dblarray*[Nthreads];
    dblarray *subDico=new dblarray[Nthreads];
    dblarray **subCode=new dblarray*[Nthreads];
    dblarray **approx=new dblarray*[Nthreads];
    dblarray *residual;
    dblarray *code;
    vector<int>* indx=new vector<int>[Nthreads]; //support

    std::cout<<"NTHREADS="<<Nthreads<<" "<<getpagesize()<<std::endl;
    code = new dblarray[Nthreads];
    residual = new dblarray[Nthreads];
    //note: getpagesize() is in bytes but I used it as nels to avoid any interference
    for(size_t kt=0;kt<Nthreads;++kt) {
        code[kt].alloc(m_natoms+getpagesize());
        residual[kt].alloc(m_npix+getpagesize());
    }
    for(size_t kt=0;kt<Nthreads;++kt){
        code[kt].reform(m_natoms);
        residual[kt].reform(m_npix);
    }

    //note: getpagesize() is in bytes but I used it as nels to avoid any interference
    for(size_t kt=0;kt<Nthreads;++kt) indx[kt].reserve(m_natoms+getpagesize());

    inputVectors=InputData->getManifoldData(m_channel);


    DEBUGME_( fits_write_dblarr("OMPDico.fits", *m_weightedDictionary); );
    //m_dictionary->display(10,27);
    //inputVectors->display(100,10);
    std::cout<<"Input Vectors="<<inputVectors->nx()<<"x"<<inputVectors->ny()<<std::endl;
    DEBUGME_( fits_write_dblarr("OMPInputVector.fits", *inputVectors); );

    #pragma omp parallel for \
        shared(inputVectors,std::cout,InputData,scalePtr,indx,\
               inputMetricPtr,allCodes,code,residual,pinvsubDico,tempname,\
                                                       subDico,subCode,approx)\
        default(none) num_threads(Nthreads) schedule(static)\
        private(stopPursuit,vectorPtr,residualPtr,residualL2Norm,new_pos,\
            codePtr,iterPtr)  reduction(+ : totalSparsity, processedVec)
    for (size_t vec=0;vec<m_ninput;vec++) {//m_ninput
        //std::cout<<" sparsecoding vector #"<<vec<<"/"<<m_ninput-1<<std::endl;
        int tid = omp_get_thread_num();

        stopPursuit=false;
        if ((m_settings.verbose == true)&&(processedVec>0))
            if(((m_ninput>100)&&((vec+1)%(m_ninput/100)==1))||(m_ninput<100))
                std::cout<<"Processing sample "<<1+vec<<" / "<< m_ninput<<
                            ", average sparsity "<<totalSparsity/processedVec<<std::endl;
        DEBUGME_( printf("PROCESS PATCH %zu\n",vec); );

        // Initialize residual = vector, index=[] and norm(residual)
        vectorPtr =InputData->getSampleChannel(vec,m_channel);
        //inputVectors->buffer()+vec*m_npix;
        residualPtr = residual[tid].buffer();
        for (size_t k=0;k<m_npix;++k,++residualPtr,++vectorPtr)
                                                     *residualPtr = *vectorPtr;
        DEBUGME_(
            if(vec==0) {
                #pragma omp critical
                {
                    sprintf(tempname[tid],"OMPresptr_%zu_%d.fits",vec,tid);
                    std::cout<<"WRITE TEMPNAME="<<tempname[tid]<<std::endl;
                    fits_write_dblarr(tempname[tid],residual[tid]);
                }
            }
        );
        residualL2Norm = OMP::inputNorm(residual[tid]); // == 0
        indx[tid].clear();

        for(size_t j=0;(j<MIN(m_settings.sparsityTarget,m_natoms))&&
                        (residualL2Norm>(1+1e-5)*m_settings.errorTarget);j++) {
            //std::cout << "   Adding atom # " << j << std::endl;
            // Initialize code
            code[tid].init(0);
            codePtr = code[tid].buffer();
            // Find maximum correlating atom
            new_pos=findClosestAtom(residual[tid].buffer(),codePtr,
                                                      inputMetricPtr,scalePtr);
            if(new_pos<0) stopPursuit=true;
            DEBUGME_(
                if(vec==0) {
                    #pragma omp critical
                    {
                        sprintf(tempname[tid],"OMPcode_%zu_%zu_%d.fits",j,vec,tid);
                        fits_write_dblarr(tempname[tid],code[tid]);
                    }
                }
            );
            //Check if this gives a previous index
            for(size_t k=0; k<indx[tid].size();++k){
                if(indx[tid][k]== new_pos) {
                    stopPursuit=true;
                    DEBUGME_(
                    std::cout<<"Atom "<<new_pos<<" appears twice. "<<
                              "Stopping the pursuit."<<std::endl;
                    );
                    break;
                }
            }
            if(residualL2Norm< this->m_settings.epsilon){
                DEBUGME_(
                    std::cout<<"Residual less than "<<
                            this->m_settings.epsilon<<". Stopping the pursuit."
                                                                    << std::endl;
                );
                stopPursuit=true;
            }
            if(stopPursuit){
                DEBUGME_(
                    std::cout << j << "atoms selected." << std::endl;
                );
                break;
            } else {
                indx[tid].push_back(new_pos); // storing new index in indx
                DEBUGME_( std::cout <<"Atom "<<new_pos<<" Selected for sample "
                                                          << vec<<std::endl; );
                subDico[tid].alloc(indx[tid].size(),m_npix,0);
                //Fill sub dictionary
                //subDico is dblarray(Nindex, Npix)
                iterPtr = subDico[tid].buffer();
                double *dictPtr;
                for(size_t kpix=0;kpix<m_npix;++kpix) {
                    dictPtr = m_weightedDictionary->buffer()+kpix* m_natoms;
                    for (size_t m=0;m<indx[tid].size();++m,++iterPtr)
                        *iterPtr= *(dictPtr +indx[tid][m]);
                }
                DEBUGME_(
                    if(vec==0){
                        #pragma omp critical
                        {
                            sprintf(tempname[tid],"OMPsubdico_%zu_%zu_%d.fits",j,vec,tid);
                            fits_write_dblarr(tempname[tid],subDico[tid]);
                        }
                    }
                );

                //subDico(n,m) = TabDico(n,indx[m]);
                // if (subDico.nx()>subDico.ny()) minT = subDico.nx()*eps;
                // else minT = subDico.ny()*eps;
                //Invert the sub-dictionary
                //tolerance estimate should be set before
                //Go_Time("computeCode::vec::pinvert");
                pinvsubDico[tid]= m_matdecomp->pinvert(subDico[tid]);
                DEBUGME_(
                    if(1){
                        #pragma omp critical
                        {
                            std::cout<<"VEC:"<<vec<<" J:"<<j<<" ASSIGN PINVSUB "<<tid
                                                    <<","<<pinvsubDico[tid]<<std::endl;
                            sprintf(tempname[tid],"OMPpinvsubdico_%zu_%zu_%d.fits",j,vec,tid);
                            fits_write_dblarr(tempname[tid],*(pinvsubDico[tid]));
                        }
                    }
                );

                //Og_Time("computeCode::vec::pinvert");
                //Compute the code
                // Computing coefficient approximation
                if(j>0) delete subCode[tid];
                subCode[tid]=m_matdecomp->matRowVec(*(pinvsubDico[tid]), *inputVectors,vec);
                DEBUGME_(
                    if(1){
                        #pragma omp critical
                        {
                            sprintf(tempname[tid],"OMPorthocode_%zu_%zu_%d.fits",j,vec,tid);
                            fits_write_dblarr(tempname[tid],*(subCode[tid]));
                        }
                    }
                );
                delete (pinvsubDico[tid]);
                // Computing vector approximation
                approx[tid]=m_matdecomp->matRowVec(subDico[tid], *(subCode[tid]),0);
                subDico[tid].free();
                vectorPtr = inputVectors->buffer()+vec*m_npix;
                iterPtr = approx[tid]->buffer();
                residualPtr=residual[tid].buffer();
                for (size_t k=0;k<m_npix;++k,++residualPtr,++vectorPtr,++iterPtr)
                    *residualPtr = *vectorPtr-*iterPtr;
                DEBUGME_(
                    if(vec==0){
                        #pragma omp critical
                        {
                            sprintf(tempname[tid],"OMPapprox_%zu_%zu_%d.fits",j,vec,tid);
                            fits_write_dblarr(tempname[tid],*(approx[tid]));
                            sprintf(tempname[tid],"OMPResUpdate_%zu_%zu_%d.fits",j,vec,tid);
                            fits_write_dblarr(tempname[tid],residual[tid]);
                        }
                    }
                );
                delete approx[tid];
                // Updating residual energy
                residualL2Norm= OMP::inputNorm(residual[tid]);
                DEBUGME_(
                if(vec<10) printf("RES NORM[%zu]=%lf\n",
                                                         vec,residualL2Norm););
            }
        }
        // Filling indx indices with coefficients
        double* allCodePtr = allCodes->buffer()+vec*m_natoms;
	    if(indx[tid].size()>0){
            codePtr = subCode[tid]->buffer();
            for (size_t k=0;k < indx[tid].size();++k,++codePtr) {
                *(allCodePtr+indx[tid][k]) = *codePtr;
                if(fabs(*codePtr)> m_settings.epsilon) totalSparsity++;
            }
            delete (subCode[tid]);
        }
        processedVec++;
    } // for vec = input vector

    delete[] indx;
    delete[] pinvsubDico;
    delete[] subDico;
    delete[] subCode;
    delete[] code;
    delete[] approx;
    delete[] residual;
    delete inputVectors;
    if (m_settings.verbose == true) std::cout << "OMP processing completed" << std::endl;
    // Time_Display_All;
    return allCodes;
}

dblarray* OMP::computeApprox(dblarray* Code, DataStruc* Input){
    return m_matdecomp-> matMultTransp(*Code,*m_dictionary);
}

void OMP::computeApprox(dblarray* Code, dblarray* Approx, DataStruc* Input){
     m_matdecomp-> matMultTransp(*Approx, *Code,*m_dictionary);
}

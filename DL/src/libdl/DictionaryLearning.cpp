/******************************************************************************
 **                   Copyright (C) 2016 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Florent Sureau [MOD rewritten from S. Beckouche Code]
 **
 **    Date:  06-07/2016
 **
 **    File:  DictionaryLearning.cpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Implementation of dictionary learning techniques.
 **    -----------
 **  - 06-07/2016: MOD, kSVD implemented
 ******************************************************************************/

#include "DictionaryLearning.h"
#include "IM_IO.h"
#include "ISAP_verbose.hpp"
#include "myexit.hpp"

#if ENABLE_VERBOSE_DEBUG_IN_DictionaryLearning
    #define DEBUGME_(x) x
#else
    #define DEBUGME_(x)
#endif

/******************************************************************************/
// SparseDL
/******************************************************************************/
SparseDL::~SparseDL() {
    DEBUGME_( printf("DELETE DICTIONARY\n"); );
    if (m_dictionary!=NULL) delete m_dictionary;
    m_dictionary=NULL;
    if (m_codetraining!=NULL) {
        DEBUGME_( printf("DELETE CODE\n"); );
        delete m_codetraining;
    }
    m_codetraining =NULL;
    if (m_approx!=NULL) {
        DEBUGME_( printf("DELETE APPROX\n"); );
        printf("DELETE APPROX\n");
        delete m_approx;
    }
    m_approx =NULL;
    if(m_dictionaryInit!=NULL) delete m_dictionaryInit;
    m_dictionaryInit = NULL;
    if(m_patchMean!=NULL) delete m_patchMean;
    m_patchMean=NULL;
}

/*!
 * Default initialization of parameters.
 */
void SparseDL::init() {
    std::cout<<"INIT SPARSEDL"<<std::endl;
    m_coder=NULL;
    m_training=NULL;
    m_dictionary=NULL;
    m_approx =NULL;
    m_codetraining =NULL;
    m_patchMean=NULL;
    m_niter=0;
    m_natoms=0;
    m_npix=0;
    m_ntrain=0;
    m_verbose=true;
    m_eps=EpsDoublePrec;
    m_atom_thresh=1;
    m_channel=0;
}

/**
 * Set Initial Dictionary.
  * @param[in] initDict Initial Dictionary to start learning with.
*/
void SparseDL::setDictionary(dblarray const &initDict) {
    FINFO_;
    if(m_dictionary!=NULL)
        delete m_dictionary;
    m_dictionary = new dblarray(initDict);
    m_coder->setDictionary(*m_dictionary);
    //! @remark AW : don't understand why this makes a copy
    //! @remark FCS: 'new' will use the constructor by copy defined in TempArray.
    FOFNI_;
}

int SparseDL::setTraining(DataStruc &training,int Channel, bool online) {
    FINFO_;
    if(online){
        std::cout<<"Online Version not implemented for this algorithm."<<
                                                "Use static instead"<<std::endl;
    }
    m_training = &training;
    if(!training.isExtracted()) training.extractPatches();
    m_ntrain = training.getNPatches();
    m_npix=training.getManifoldDimSize(Channel);
    DEBUGME_(std::cout << " setTraining : " << m_ntrain << " atoms of "
                                          << m_npix << " pixels" << std::endl;);
    if(m_patchCentering)
        centerPatches();
    if(m_trainingSetCentering)
        centerTraining();
    FOFNI_;
    return 0;
}

//! Create a random initial dictionary
int SparseDL::initDictionary(size_t nAtoms, double density,
                                            InitDictionaryType dictionaryType) {
    FINFO_;
    // Create Initial Dictionary
    std::cout<<"INNER m_natoms="<<nAtoms<<" CHAN="<<m_channel<<std::endl;
    m_natoms = nAtoms;
    if(m_training == NULL) {
        std::cout << "Error: initDictionary can only be called"<<
                                             " after setTraining" << std::endl;
        return -1;
    }
    switch(dictionaryType) {
        case InitDictionaryType::UserDefinedDictionary:
            if (m_verbose) std::cout << "Using a user-defined "<<
                                             "initial dictionary. " <<std::endl;
            return 0;
            break;
        case InitDictionaryType::RandomAtomsDictionary:
            if (m_verbose) std::cout << "Set initial dictionary from random training samples. " <<std::endl;
            dblarray* TempData;
            TempData=m_training->getManifoldData(m_channel);
            m_dictionaryInit = new RandomDictInit(nAtoms, TempData);
            delete TempData;
            break;
        case InitDictionaryType::RandomSparseDictionary:
            if (m_verbose) std::cout << "Set initial dictionary as random sparse with density: " << density <<std::endl;
            m_dictionaryInit = new SparseDictInit(density, nAtoms, m_npix);
            break;
    }
    m_dictionaryInit->setRNG(m_RNG);
    dblarray *Dict = m_dictionaryInit->computeInit();
    setDictionary(*Dict);
    delete Dict;
    FOFNI_;
    return 0;
}

void SparseDL::setPatchCentering(bool patchCentering) {
    if(m_training != NULL) {
        std::cerr << "SparseVectorDL::setPatchCentering can only be called before initializing the training set" << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    else
        m_patchCentering = patchCentering;
}
void SparseDL::setTrainingSetCentering(bool trainingSetCentering) {
    if(m_training != NULL) {
        std::cerr << "SparseVectorDL::setTrainingSetCentering can only be called before initializing the training set" << std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    else
        m_trainingSetCentering = trainingSetCentering;
}

/******************************************************************************/
// SparseVectorDL
/******************************************************************************/


//Setters
void SparseVectorDL::setDictionary(dblarray &initDict){
    SparseDL::setDictionary(initDict);
    m_natoms=m_dictionary->nx();
    m_npix=m_dictionary->ny();
}



int SparseVectorDL::setTraining(DataStruc &training,int Channel, bool online) {
    FINFO_;
    SparseDL::setTraining(training,Channel,online);
    setChannel(Channel);
    FOFNI_;
    return 0;
}
//Generic Routines

//Each patch is centered
void SparseVectorDL::centerPatches(){
    FINFO_;
    if(getManifold()==Manifold::Rn){
        std::cout <<"CENTER THE PATCHES"<<std::endl;
        dblarray* tempData=m_training->getManifoldData(m_channel);
        if(m_patchMean!=NULL) delete m_patchMean;
        this->m_matdecomp->centerRow2DData(*tempData,m_patchMean);
        m_training->putManifoldData(tempData,m_channel);
        delete tempData;
        m_patchCentering =true;
    } else if(getManifold()==Manifold::Rpn)
        std::cout<<"DO NOT PERFORM centerPatches for Rpn"<<std::endl;
    else {
        std::cout<<"NOT RECOGNIZED centerPatches manifold"<<std::endl;
        exit(EXIT_FAILURE);
    }
    FOFNI_;
}


//Training set is centered and mean patch saved
void SparseVectorDL::centerTraining(char *nameMeanPatch){
    FINFO_;
    if(getManifold()==Manifold::Rn){
        dblarray* tempData=m_training->getManifoldData(m_channel);
        this->m_matdecomp->centerCol2DData(*tempData, nameMeanPatch);
        m_training->putManifoldData(tempData,m_channel);
        delete tempData;
    } else if((getManifold()==Manifold::S1n)||(getManifold()==Manifold::S1Rn)
                                            ||(getManifold()==Manifold::S1wRn)){
        std::cout<<"CENTER TRAINING FOR S1n NOT YET IMPLEMENTED"<<std::endl;
        MyExit::exit(EXIT_FAILURE);
    } else {
        std::cout<<"CENTER TRAINING NOT YET IMPLEMENTED"<<std::endl;
        MyExit::exit(EXIT_FAILURE);

    }
    FOFNI_;
}

//Center dictionary
void SparseVectorDL::centerDictionary(){
    FINFO_;
    this->m_matdecomp->centerCol2DData(*m_dictionary);
    FOFNI_;
}

//Project training patch into l2 ball
void SparseVectorDL::l2projTraining(){
    dblarray* tempData=m_training->getManifoldData(m_channel);
    this->m_matdecomp->l2ProjRow2DData(*tempData);
    m_training->putManifoldData(tempData,m_channel);
    delete tempData;
}

//Project each atom into l2 ball
void SparseVectorDL::l2projDictionary(){
    this->m_matdecomp->l2ProjCol2DData(*m_dictionary);
}

//Project and center each training patch into l2 ball
void SparseVectorDL::normTraining(){
    dblarray* tempData=m_training->getManifoldData(m_channel);
    this->m_matdecomp->l2NormalizeRow2DData(*tempData);
    m_training->putManifoldData(tempData,m_channel);
    delete tempData;
}

//Project and center each atom into l2 ball
void SparseVectorDL::normDictionary(){
    this->m_matdecomp->l2NormalizeCol2DData(*m_dictionary);
}

//Replace unused or not so used atoms (m_atom_thresh) by random training samples
unsigned int SparseVectorDL::replaceUnusedAtoms(){
    size_t atom_usage,atoms_replaced, new_sample_ind;
    double atom_usageNorm;
    dblarray atom(m_npix);
    double *atomBuf,*tempBuf,*iterBuf;
    //Replace unused atoms by selecting randomly a training example
    //(or a training sample badly represented for S1n)
    atoms_replaced = 0;
    tempBuf=m_codetraining->buffer();
    DEBUGME_( m_RNG->setSeed(19800706); );
    std::list<unsigned int> listNoDecomp;
    std::list<unsigned int>::iterator itList;
    size_t ka;
    for(size_t kv=0;kv<m_ntrain;++kv){
        for(ka=0;ka<m_natoms;++ka)
            if(fabs((*m_codetraining)(ka,kv))>0) break;
        if(ka==m_natoms) listNoDecomp.push_back(kv);
    }
    std::cout<<"NON CODED SAMPLES="<<listNoDecomp.size();

    for(size_t na=0;na<m_natoms;++na,++tempBuf){
        atom_usage = 0;
        atom_usageNorm = 0.f;
        iterBuf=tempBuf;
        for (size_t p=0;p<m_ntrain;++p, iterBuf+=m_natoms)
            if(fabs(*iterBuf)>m_eps) atom_usage++;
        //std::cout<<"Atom["<<na<<"] usage "<<atom_usage<<" - "<<m_eps<<std::endl;
        if (atom_usage < m_atom_thresh){ //Atom is not sufficiently used
            if((getManifold()==Manifold::S1n)||(getManifold()==Manifold::S1wRn)){
                if(!listNoDecomp.empty()){
                    new_sample_ind= m_RNG->getUniformInt(listNoDecomp.size());
                    itList = listNoDecomp.begin();
                    std::advance(itList, new_sample_ind);
                    new_sample_ind= *itList;
                    itList = listNoDecomp.erase (itList);
                } else new_sample_ind= m_RNG->getUniformInt(m_ntrain);//prev m_ntrain-1?
            } else new_sample_ind= m_RNG->getUniformInt(m_ntrain);//prev m_ntrain-1?

            DEBUGME_( printf("Replace atom=%zu used %zu times by %zu \n",
                                             na, atom_usage, new_sample_ind); );
            atomBuf=atom.buffer();
            getTrainingSample(new_sample_ind, atomBuf);
            if(getManifold()==Manifold::Rn){
                if(m_patchCentering) this->m_matdecomp->l2NormalizeRow2DData(atom);
                else this->m_matdecomp->l2ProjRow2DData(atom);
            } else if (getManifold()==Manifold::Rpn){
                this->m_matdecomp->l2ProjRow2DData(atom);
            }
            atomBuf=atom.buffer();
            this->setAtom(na, atomBuf,1);
            atoms_replaced++;
        } else {
            DEBUGME_( printf("Keeping atom=%zu used %zu times, norm=%f,mean=%f\n", na, atom_usage,atom_usageNorm,atom_usageNorm/atom_usage); );
        }
    }
    std::cout<<"NEW NON CODED SAMPLES:"<<listNoDecomp.size()<<std::endl;
    return atoms_replaced;
}

//Replace non significant or highly correlated (=) atoms
unsigned int SparseVectorDL::cleanDictionary(){
    DEBUGME_( char tempname[1024]; );
    size_t atoms_replaced=0, new_sample_ind, n_loop =0;
    //If manifold is S1n, then we do not want to normalize atoms (they are
    //already normalized in S1n (or in [-pi,pi[ if we have an intrinsic
    //representation))
    if(getManifold()==Manifold::Rn) this->l2projDictionary();
    else if(getManifold()==Manifold::Rpn) this->l2projDictionary();
    dblarray *gram= this->m_matdecomp->gram(*m_dictionary);
    gsl_permutation* Perm;
    dblarray atom(m_npix,1),*scProd;
    double *gramBuf=gram->buffer(),*atomBuf,*bufSc,maxIP,minl2Diff,l2Diff;
    double currEnergy;
    size_t kdiff;
    DEBUGME_({ sprintf(tempname,"DLDico_inputgram.fits");
               fits_write_dblarr(tempname,*m_dictionary);
               sprintf(tempname,"DLGramDico.fits");
               fits_write_dblarr(tempname,*gram);
             });
    atoms_replaced = 0;

    std::vector<unsigned int> listNoDecomp;
    if(((getManifold()==Manifold::S1n)||(getManifold()==Manifold::S1wRn))
                                                    &&(m_codetraining!=NULL)){
        size_t ka;
        for(size_t kv=0;kv<m_ntrain;++kv){
            for(ka=0;ka<m_natoms;++ka)
                if(fabs((*m_codetraining)(ka,kv))>0) break;
            if(ka==m_natoms) listNoDecomp.push_back(kv);
        }
    }
    size_t na;
    for(na=0;na<m_natoms;++na){
        n_loop=0;
        Perm= m_RNG->createPermutation(m_ntrain);//prev m_ntrain-1?
        atomBuf=atom.buffer();
        getAtom(na, atomBuf);
        while(n_loop<m_ntrain) {
            //std::cout << "GRAM["<<na<<","<<na<<"]="<<gramBuf[na+na*m_natoms]<<std::endl;
            if((getManifold()==Manifold::Rn)||(getManifold()==Manifold::Rpn)){
                //Replace highly correlated atoms and 0 atoms
                if(gramBuf[na+na*m_natoms]>1e-6){
                    scProd=this->m_matdecomp->matMult(atom,*m_dictionary);
                    bufSc=scProd->buffer();
                    maxIP=0.;
                    kdiff=0;
                    for(size_t ka=0;ka<na;++ka,++bufSc){
                        gramBuf[na+ka*m_natoms]=* bufSc;
                        if(fabs(*bufSc)>maxIP){
                            maxIP=fabs(*bufSc);
                            kdiff=ka;
                        }
                    }
                    for(size_t ka=na+1;ka<m_natoms;++ka){
                        ++bufSc;
                        gramBuf[ka+na*m_natoms]=*bufSc;
                        if(fabs(*bufSc)>maxIP){
                            maxIP=fabs(*bufSc);
                            kdiff=ka;
                        }
                    }
                    delete scProd;
                    if(maxIP<0.999999) {
                        setAtom(na, atom.buffer(),1);
                        break;
                    }
                }
            } else {
                //We accept 0 atoms, we accepted unnormalized atoms, we look
                //instead at the closeness of two atoms
                scProd=this->m_matdecomp->matMult(atom,*m_dictionary);
                //\theta_i\theta_j
                bufSc=scProd->buffer();
                currEnergy=atom.energy();
                minl2Diff =currEnergy+MAX((*gram)(0,0),(*gram)(1,1));
                kdiff=0;
                //check for first atoms what is the l2 diff
                for(size_t ka=0;ka<na;++ka,++bufSc){
                    l2Diff= currEnergy+gramBuf[ka+ka*m_natoms]-2*(*bufSc);
                    if(l2Diff<0.) std::cout<<"Problem l2 Diff:"<<l2Diff<<std::endl;
                    //theta_i^2+theta_j^2-2*theta_i*theta_j
                    //std::cout << "GRAM["<<na<<","<<ka<<"]="<<*bufSc<<std::endl;
                    if(l2Diff<minl2Diff){
                        minl2Diff=l2Diff;
                        kdiff=ka;
                    }
                    gramBuf[na+ka*m_natoms]=* bufSc;
                }
                for(size_t ka=na+1;ka<m_natoms;++ka){
                    ++bufSc;
                    l2Diff= currEnergy+gramBuf[ka+ka*m_natoms]-2*(*bufSc);
                    if(l2Diff<0.) std::cout<<"Problem l2 Diff:"<<l2Diff<<std::endl;
                    //std::cout << "GRAM["<<na<<","<<ka<<"]="<<*bufSc<<std::endl;
                    if(l2Diff<minl2Diff){
                        minl2Diff=l2Diff;
                        kdiff=ka;
                    }
                    gramBuf[ka+na*m_natoms]=*bufSc;
                }
                delete scProd;
                if(minl2Diff > 0.000001) {//before: 0.000001
                    setAtom(na, atom.buffer(),1);
                    break;
                }
            }
            if((!listNoDecomp.empty())&&(listNoDecomp.size()> n_loop)){
                new_sample_ind= listNoDecomp[n_loop];
            } else {
                std::cout<<"Atom "<<na<< ": matching "<<kdiff<<"("<<minl2Diff<<")"<<std::endl;
                new_sample_ind =gsl_permutation_get(Perm,n_loop);
            }
            std::cout<<"Atom "<<na<< ": new_sample_ind "<<new_sample_ind <<std::endl;
            atomBuf=atom.buffer();
            getTrainingSample(new_sample_ind, atomBuf);
            if(getManifold()==Manifold::Rn){
                if(m_patchCentering) this->m_matdecomp->l2NormalizeRow2DData(atom);
                else this->m_matdecomp->l2ProjRow2DData(atom);
            } else if(getManifold()==Manifold::Rpn)
                                    this->m_matdecomp->l2ProjRow2DData(atom);

            gramBuf[na+na*m_natoms]=atom.energy();
            ++n_loop;
        }
        if(n_loop>0) ++atoms_replaced;
        if(n_loop==m_ntrain){
            char tempnameFail[1024];
            std::cout << "Cannot find appropriate training sample" << std::endl;
            sprintf(tempnameFail,"DLDico_cannotFindNewSample.fits");
            fits_write_dblarr(tempnameFail,*m_dictionary);
            MyExit::exit(EXIT_FAILURE);
        }
    }
    delete gram;
    return atoms_replaced;
}

/******************************************************************************/
//DL using full batch training set
/******************************************************************************/
//Initialization of DL
int SparseVectorDL:: initDL(){
    DEBUGME_( char tempname[1024]; );
    unsigned int replaced_atoms;
    DEBUGME_({ sprintf(tempname,"DLDico.fits");
               fits_write_dblarr(tempname,*m_dictionary);
             });
    DEBUGME_({ sprintf(tempname,"DLTraining.fits");
               fits_write_dblarr(tempname,* m_training->getDataPtr());
             });
    if(this->m_coder==NULL){
        std::cout << "SparseVectorDL::Coder not initialized" <<std::endl;
        return -1;
    }
    if(this->m_dictionary ==NULL){
        std::cout << "SparseVectorDL::Dictionary not initialized" <<std::endl;
        return -1;
    }
    if(((size_t) this->m_dictionary->nx()!=m_natoms)||
       ((size_t) this->m_dictionary->ny()!=m_npix)){
        std::cout << "SparseVectorDL::initDL : Wrong size dictionary: [" <<
                     this->m_dictionary->nx() <<","<<this->m_dictionary->ny()<<"] vs ["
                  <<m_natoms<<","<<m_npix<<"]"<<std::endl;
        return -1;
    }
    if(this->m_training ==NULL){
        std::cout << "SparseVectorDL::Training set not initialized" <<std::endl;
        return -1;
    }
    if(((size_t) this->m_training->getManifoldDimSize(m_channel)!=m_npix)||
       ((size_t) this->m_training->getNPatches() != m_ntrain)){
        std::cout << "SparseVectorDL::Wrong size training set: [" <<
                     this->m_training->getManifoldDimSize(m_channel) <<","<<
                     this-> m_training->getNPatches()<<"] vs ["
                     <<m_npix<<","<<m_ntrain<<"]"<<std::endl;
        return -1;
    }
    if(m_approx!=NULL) delete m_approx;
    m_approx=NULL;
    //l2 proj dictionary
    if(getManifold()==Manifold::Rn) {
        if(m_patchCentering) this->normDictionary();
        else this->l2projDictionary();
    } else {
        std::cout<<"INIT DL: manifold not recognized"<<std::endl;
        exit(EXIT_FAILURE);
    }
    replaced_atoms=cleanDictionary();
    //set a reference to the dictionary for the coder
    this->m_coder->setDictionary(*m_dictionary); //! @remark AW : this seems to be a copy, not a pointer. bug ?
    //! @remark FCS You're right, it should be a pointer here as the sparse
    // coder is not supposed to modify the dictionary.
    if(m_verbose)
        std::cout <<"SparseVectorDL::Replaced atoms for init Dictionary: " <<
                    replaced_atoms << std::endl;

    DEBUGME_({ sprintf(tempname,"DLDicoProjSPAM.fits");
               fits_write_dblarr(tempname,*m_dictionary);
             });

    if (m_verbose) {
        this->m_coder->setVerbose(true);
        std::cout << "SparseVectorDL::Learning dictionary of" << m_natoms <<
                     " atoms of " << m_npix <<" pixels" << std::endl;
        std::cout << "Using DL " << this->getSignature ()<< std::endl;
        std::cout <<"SparseVectorDL::Using Matrix Decomposition: " <<
                    m_matdecomp->getSignature() << std::endl;
        std::cout << "SparseVectorDL::Starting Dictionary Learning for " <<
                     m_niter  << " iterations" << std::endl;
        std::cout << "Using coder " << m_coder->getSignature() <<
                     " with SparsityTarget = "<< m_coder-> getSparseTarget() <<
                     " and ErrorTarget = " << m_coder->getErrorTarget() << endl;
    }
    return 0;
}

double SparseVectorDL::computeAverageError(bool recompute){
    double average_error = 0;
    double *sparseBuf, *trainBuf;
    if(m_approx !=NULL) delete m_approx;
    m_approx=getApprox(recompute);
    sparseBuf= m_approx->buffer();
    if(getManifold()==Manifold::Rn){
        for(size_t kpatch=0;kpatch<m_ntrain;++kpatch){
            trainBuf=m_training->getSampleChannel(kpatch,m_channel);
            for (size_t kel=0;kel<m_npix;++kel,++sparseBuf,++trainBuf)
                average_error += (*sparseBuf - *trainBuf)*(*sparseBuf - *trainBuf);
        }
    }  else {
        std::cout<<"CANNOT COMPUTE AVERAGE ERROR. Manif not recognized"<<std::endl;
        exit(EXIT_FAILURE);
    }
    average_error = sqrt(average_error)/m_ntrain;
    return average_error;
}

//Main learning procedure for vectorized data
int SparseVectorDL::learning(){
    //Average results
    double average_sparsity = 0;
    double average_error = 0;
    //Code
    //atoms
    unsigned int atoms_replaced;
    DEBUGME_( char tempname[1024]; );

    if(initDL()==-1){
        std::cout << "DL Initialization failed" << endl;
        MyExit::exit(EXIT_FAILURE);
    }
    if(getManifold()!=m_coder->getManifold()){
        std::cout<<"BEWARE: INCONSISTENCY BETWEEN DL MANIFOLD "<<getManifold()
                <<" AND CODER MANIFOLD "<< m_coder->getManifold()<<std::endl;
    }
    if(getManifold()!=m_training->getManifold(m_channel)){
        std::cout<<"BEWARE: INCONSISTENCY BETWEEN DL MANIFOLD "<<getManifold()
                                <<" AND TRAINING SET MANIFOLD "
                                <<m_training->getManifold(m_channel)<<std::endl;
    }
    for (size_t it=0;it<m_niter;++it){
        //Print some informative stuff
        if (m_verbose){
            if (m_niter > 100){
                if ((it +1)%(m_niter/10) == 1)
                    std::cout << "DL iteration " << it +1 << " / " << m_niter
                              << ", average sparsity " << average_sparsity <<
                                 ", average error " << average_error << std::endl;
            }
            else std::cout << "DL iteration " << it +1 << " / " << m_niter <<
                              ", average sparsity " << average_sparsity <<
                              ", average error " << average_error << std::endl;
        }
        //Sparse Coding
        if(m_codetraining!=NULL) delete m_codetraining;
        m_codetraining= m_coder->computeCode(m_training);//
        DEBUGME_({ sprintf(tempname,"DLCodePrior_it%zu.fits",it);
                   fits_write_dblarr(tempname,*m_codetraining);
                 });
        // Computing average sparsity given sparse encoding coefficients
        if(m_verbose){
            average_sparsity = 0;
            double *TrainingBuffer= m_codetraining->buffer();
            for(size_t kel=0;kel<m_ntrain*m_natoms;++kel,++TrainingBuffer)
                if(*TrainingBuffer !=0.0) average_sparsity++;
            average_sparsity = average_sparsity / (m_ntrain);
            std::cout << "Sparse coding completed, average sparsity "
                      << average_sparsity <<std::endl;
            std::cout << "Updating dictionary ..." << std::endl;
            average_error=computeAverageError(false);//No need to recompute the code
            DEBUGME_(if(1){ std::cout << "Average Error prior update:"<< average_error <<std::endl;
                            sprintf(tempname,"DLDictPrior_it%zu.fits",it);
                            fits_write_dblarr(tempname,*m_dictionary);
                            sprintf(tempname,"DLApproxPrior_it%zu.fits",it);
                            fits_write_dblarr(tempname,*m_approx);
                     });
            delete m_approx;
            m_approx=NULL;
        }
        //MOD Dictionary Update
        this->updateDictionary();
        // Computing sparse approximation and average quadratic error
        if (m_verbose){
            average_error=computeAverageError(true);//Need to recompute the code
            DEBUGME_(if(1){  std::cout << "Average Error post update:"<< average_error <<std::endl;
                            sprintf(tempname,"DLUpDict_it%zu.fits",it);
                            fits_write_dblarr(tempname,*m_dictionary);
                            sprintf(tempname,"DLUpApprox_it%zu.fits",it);
                            fits_write_dblarr(tempname,*m_approx);
                            sprintf(tempname,"DLUpCode_it%zu.fits",it);
                            fits_write_dblarr(tempname,*m_codetraining);
                     });
        }
        atoms_replaced = replaceUnusedAtoms();
        atoms_replaced += cleanDictionary();
        delete m_codetraining;
        delete m_approx;
        m_approx=NULL;
        m_codetraining=NULL;
        if ((m_verbose) && (atoms_replaced !=0)){
            if (m_niter >100) {
                if ((it+1)%((int)m_niter/10) == 1)
                    std::cout << "Replaced " << atoms_replaced <<
                                 " unused atoms with random training samples" << std::endl;
            }  else
                std::cout << "Replaced " << atoms_replaced <<
                             " unused atoms with random training samples" << std::endl;
        }
        //If manifold is S1n, then we do not want to normalize atoms (they are
        //already normalized in S1n (or in [-pi,pi[ if we have an intrinsic
        //representation))
        if((getManifold()==Manifold::Rn)||(getManifold()==Manifold::Rpn))
            this->m_matdecomp->l2ProjCol2DData(*m_dictionary);
        else {
            std::cout<<"LEARNING : Manifold not recognized"<<std::endl;
            exit(EXIT_FAILURE);
        }
        // Updating sparse coder with new version of dictionary
        m_coder->setDictionary(*m_dictionary); //! @remark AW why not use a pointer ?

    }
    return 0;
}

/******************************************************************************/
//MOD
/******************************************************************************/
//Implementation of dictionary update for MOD for vectorized data
void SparseVectorMOD::updateDictionary(){
    dblarray* codePinv;
    double minThresh;
    DEBUGME_( char tempname[1024]; );
    minThresh = (double) m_ntrain* m_eps* sqrt(m_codetraining->energy());
    DEBUGME_( printf("Minthresh=%lf \n", minThresh); );
    //Chosing eigenvalue threshold value
    if(getManifold()==Manifold::Rn){
        if(m_dictionary!=NULL) delete m_dictionary;
        codePinv=m_matdecomp->pinvert(*m_codetraining,minThresh);
        DEBUGME_({ printf("CodePinv, Nx=%d cols Ny=%d rows\n", codePinv->nx(), codePinv->ny());
                   sprintf(tempname,"MODPinvCode.fits");
                   fits_write_dblarr(tempname,*codePinv);
                 });
        dblarray* TempArray=m_training->getManifoldData(m_channel);
        m_dictionary=m_matdecomp->transpMatMult(*codePinv,*TempArray);
        delete codePinv;
        delete TempArray;//Fake structure if Data Structure is contiguous
        if(getManifold()==Manifold::S1Rn){
            double *DicoBuf= m_dictionary->buffer();
            for(size_t kel=0;kel<(size_t) m_dictionary->n_elem();++kel,
                                                                    ++DicoBuf){
                if(*DicoBuf>M_PI) *DicoBuf=M_PI;
                if(*DicoBuf<-M_PI) *DicoBuf=-M_PI;
            }
        }
    }
}

/******************************************************************************/
//K-SVD and Approximate K-SVD
/******************************************************************************/
//Create the list of training samples used by a given atom and list of coefs
size_t SparseVectorKSVD::createListSubMatrix(size_t AtomNb,size_t*listTrain,
                                             double* listCoefs){
    double *iterBuf = m_codetraining->buffer()+AtomNb;

    size_t nlistTrain=0;
    for(size_t nt=0;nt<m_ntrain;++nt, iterBuf +=m_natoms){
        if(fabs(*iterBuf)>0) {
            listTrain[nlistTrain]= nt;
            listCoefs[nlistTrain]=*(iterBuf);
            ++nlistTrain;
        }
    }
    return nlistTrain;
}

//Implementation of dictionary update for KSVD/approx KSVD for vectorized data
void SparseVectorKSVD::updateDictionary(){
    size_t* listTrain=(size_t*)malloc(m_ntrain*sizeof(size_t));
    double* listCoefs=(double*)malloc(m_ntrain*sizeof(double));
    dblarray * subCodeMatrix;
    double * subCodeMatrixBuf;
    //compute approximation
    if(m_approx !=NULL) delete m_approx;
    m_approx=m_matdecomp->matMultTransp(*m_codetraining,*m_dictionary);
    //compute full residual
    dblarray residual(m_npix,m_ntrain);
    double* residualBuf= residual.buffer();
    double *trainBuf;

    double *approxBuf= m_approx->buffer();
    for(size_t kpatch=0;kpatch<m_ntrain;++kpatch){
        trainBuf=m_training->getSampleChannel(kpatch,m_channel);
        for(size_t kel=0;kel<m_npix;++kel,++residualBuf,++trainBuf,
            ++ approxBuf)
            *residualBuf= (*trainBuf)-(* approxBuf);
    }
    delete m_approx;
    m_approx=NULL;
    size_t NkeptTrain;
    dblarray Atom(m_npix);
    double *iterAtom=Atom.buffer(),*iterCoefs,*VtBuf;
    double maxEigen,AtomEnergy;
    dblarray U,S,Vt;

    for(size_t ka=0;ka<m_natoms;ka++) {
        getAtom(ka,Atom.buffer());
        //Get list of examples involving column ka
        NkeptTrain= createListSubMatrix(ka, listTrain,listCoefs);
        if(NkeptTrain>0) {
            subCodeMatrix=new dblarray(NkeptTrain,m_npix);
            subCodeMatrixBuf= subCodeMatrix->buffer();
            iterCoefs= listCoefs;
            for(size_t ks=0;ks<NkeptTrain;++ks,++iterCoefs) {
                //go to selected examples
                residualBuf= residual.buffer()+listTrain[ks]*m_npix;
                iterAtom= Atom.buffer();
                subCodeMatrixBuf=subCodeMatrix->buffer()+ks;
                for(size_t kpix=0;kpix<m_npix;++ kpix, subCodeMatrixBuf+=NkeptTrain,
                    ++residualBuf,++iterAtom){
                    (*residualBuf)+=(*iterCoefs)*(*iterAtom);//Update residual
                    *(subCodeMatrixBuf)=(*residualBuf);//Construct sub-matrix
                }
            }
            if(m_approxSVD) {
                for(size_t kit=0;kit<m_approxSVDIter;++kit){//Inner Loop
                    //Update Atom
                    Atom.init(0.);
                    AtomEnergy=0.;
                    iterAtom=Atom.buffer();
                    for(size_t kpix=0;kpix<m_npix;++ kpix,++iterAtom){
                        iterCoefs= listCoefs;
                        subCodeMatrixBuf=subCodeMatrix->buffer()+kpix* NkeptTrain;
                        for(size_t ks=0;ks<NkeptTrain;++ks,++iterCoefs,++subCodeMatrixBuf)
                            *iterAtom+=(*iterCoefs)*(*subCodeMatrixBuf);
                        AtomEnergy+=(*iterAtom)*(*iterAtom);
                    }
                    iterAtom=Atom.buffer();
                    AtomEnergy=sqrt(AtomEnergy);
                    for(size_t kpix=0;kpix<m_npix;++ kpix,++iterAtom)
                        *iterAtom/=AtomEnergy;
                    setAtom(ka, Atom.buffer(),1);
                    //Update Coefs
                    iterCoefs= listCoefs;
                    for(size_t ks=0;ks<NkeptTrain;++ks,++iterCoefs) {
                        subCodeMatrixBuf=subCodeMatrix->buffer()+ks;
                        iterAtom=Atom.buffer();
                        *(iterCoefs)=0;
                        for(size_t kpix=0;kpix<m_npix;++ kpix,++iterAtom,
                            subCodeMatrixBuf+=NkeptTrain)
                            *(iterCoefs)+= (*iterAtom)*(*subCodeMatrixBuf);
                    }
                }
                //set final coefs
                for(size_t ks=0;ks<NkeptTrain;++ks){
                    iterCoefs=m_codetraining->buffer()+ka+listTrain[ks]*m_natoms;
                    *iterCoefs= listCoefs[ks];
                }
            } else {
                //Compute SVD
                m_matdecomp->svd(*(subCodeMatrix),U,S,Vt);
                //Set Atom
                setAtom(ka, U.buffer(),U.nx());
                //Set Coefs
                maxEigen=S(0);
                VtBuf =Vt.buffer();
                for(size_t ks=0;ks<NkeptTrain;++ks,++VtBuf) {
                    iterCoefs=m_codetraining->buffer()+ka+listTrain[ks]*m_natoms;
                    *iterCoefs=(*VtBuf)*maxEigen;
                    listCoefs[ks]=*iterCoefs;
                }
            }
            getAtom(ka,Atom.buffer());
            //Update residual
            iterCoefs= listCoefs;
            for(size_t ks=0;ks<NkeptTrain;++ks,++ iterCoefs) {
                //go to selected examples
                residualBuf= residual.buffer()+listTrain[ks]*m_npix;
                iterAtom= Atom.buffer();
                for(size_t kpix=0;kpix<m_npix;++ kpix,++iterAtom,++residualBuf)
                    (*residualBuf)-=(*iterCoefs)*(*iterAtom);//subtract thex explained part
            }
            delete subCodeMatrix;
        }//End number of patches >0
    }
    free(listTrain);
    free(listCoefs);
}

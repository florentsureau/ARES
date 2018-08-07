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
 **    Date:  2016/07
 **
 **    File:  Dl_test_learning.cpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION   Main program for single-valued patch-based learning techniques
 **    -----------
 **
 **
 ******************************************************************************/
#include <cmath>
#include "Array.h"
#include "NR.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_Block.h"
#include <iomanip>
#include <sys/time.h>
#include "gslRNG.h"
#include "DictionaryInit.h"
#include "DictionaryLearning.h"

#ifdef _OPENMP
#include "OpenMPInterface.h"
#endif
#define _DEBUG_LEARNING
/****************************************************************************/
// 2D double array, containing dictionary with Npix rows and Natoms columns,
// stored in row major order.
// Beware of dblarray conventions: Index kpix,katom of Dico accessed as:
// Dico(katom,kpix)=buffer[katom+kpix*Natoms];
char nameDicoIn[1024];

// 2D double array, containing training set with Npix rows and Npatches cols ,
// stored in row major order, as provided by dl_extract_patches.
// Beware of dblarray conventions: Index kpix,ktrain of TrainSet accessed as:
// TrainSet(kpix,kpatch)=buffer[kpix+ kpatch*Npix];
// This allows to access fastly to all pixels of a given patch (e.g. for inner
// products)
char nameVectIn[1024];  // 2D array

//output dictionary, same convention as input
char nameDicoOut[1024];

//optional output approximation for patches, same convention as input patches.
char nameApproxOut[1024];
//optional mean patch from training set.
char nameMeanPatch[1024];
//optional input metric to use for sparse coding.
char nameInputMetric[1024];
//option reference for S1 log-map computation
char nameInputReferenceS1[1024];

//optional code for patches, with Npatches rows and Npix columns, stored in row
//major order.
// Beware of dblarray conventions: Index katom,ktrain of CodeSet accessed as:
// Codeset(katom,kpatch)=buffer[katom+ kpatch*Natoms];
// This allows to access fastly to all coefs of a given patch (e.g. for sparse
// coding)
char nameCodeOut[1024];

bool verbose=false;
bool timing=false;
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, const char *opts);
//double EpsCompare=1e-15;//Expected accuracy !! Not used, redefined in SparseDec

double ErrorTarget = -1.;
int NbIter = 10;
size_t NAtoms=0;
int SparsityTarget = -1;
int minimalSparsity = 0;
double minimalCorrelation = 0;
bool UseSparsity=true;
bool ReadInitDico=false;
double InitDensity=0.1;
bool SaveCode=false;
bool SaveApprox=false;
unsigned long int GeneratorSeed=0;
bool useSparseInit=false;
bool useRandomTrainInit=false;

unsigned int KSVDit=1;
bool UseKSVD=false;
bool UseApproxKSVD=false;
bool CenterPatches=false;
bool CenterTraining=false;
bool UseInputMetric =false;

/*********************************************************************/
static void usage(char *argv[])
{

    fprintf(OUTMAN, "Usage: %s options nameVectIn.fits nameDicoOut.fits \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    fprintf(OUTMAN, "      [-e ErrorTarget]\n");
    fprintf(OUTMAN, "      [-s TargetSparsity]\n");
    fprintf(OUTMAN, "      [-i Iteration number, default 10]\n");
    fprintf(OUTMAN, "      [-d InitialDictionary]\n");
    fprintf(OUTMAN, "      [-m name of fits file containing input metric]\n");
    fprintf(OUTMAN, "      [-n Natoms if no input dictionary]\n");
    fprintf(OUTMAN, "      [-p Name of mean patch to save and to subtract\n");
    fprintf(OUTMAN, "          from training set]\n");
    fprintf(OUTMAN, "      [-r Init Sparse Density if no input dictionary]\n");
    fprintf(OUTMAN, "      [-R Random training examples for initialization]\n");
    fprintf(OUTMAN, "      [-g generator seed - default 0]\n");
    fprintf(OUTMAN, "      [-c Name for training set code to save]\n");
    fprintf(OUTMAN, "      [-a Name for training set approx to save]\n");
    fprintf(OUTMAN, "      [-C center training patches (mean of each patch=0)]\n");
    fprintf(OUTMAN, "      [-K use KSVD]\n");
    fprintf(OUTMAN, "      [-k Iterations for approx KSVD (used), default 1]\n");
    fprintf(OUTMAN, "      [-M Minimal Correlation]\n");
    fprintf(OUTMAN, "      [-N Minimal Sparsity]\n");
    fprintf(OUTMAN, "      [-T Timing on]\n");
    fprintf(OUTMAN, "      [-V verbose]\n");
    exit(EXIT_SUCCESS);
}
/*********************************************************************/
/* GET COMMAND LINE ARGUMENTS */
static void transinit(int argc, char *argv[])
{
    int c,outer_loop_threads,inner_loop_threads;
    /* get options */
    while ((c = GetOpt(argc,argv,"a:c:d:e:f:g:i:k:m:n:p:r:s:t:CKMNRTV")) != -1){
        switch (c){
            case 'a':
            if (sscanf(OptArg,"%s", nameApproxOut) != 1){
                fprintf(OUTMAN, "bad -a value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            SaveApprox =true;
            break;
            case 'c':
            if (sscanf(OptArg,"%s", nameCodeOut) != 1){
                fprintf(OUTMAN, "bad -c value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            SaveCode=true;
            break;
            case 'd':
            if (sscanf(OptArg,"%s", nameDicoIn) != 1){
                fprintf(OUTMAN, "bad -d value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            ReadInitDico=true;
            break;
            case 'e':
            if(sscanf(OptArg,"%lf", &ErrorTarget) !=1){
                fprintf(OUTMAN, "bad -e value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            break;
            case 'g':
            if(sscanf(OptArg,"%lu", &GeneratorSeed) !=1){
                fprintf(OUTMAN, "bad -g value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            break;
            case 'i':
            if(sscanf(OptArg,"%d", &NbIter) !=1) {
                fprintf(OUTMAN, "bad -i value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            break;
            case 'k':
            if (sscanf(OptArg,"%u", &KSVDit) != 1){
                fprintf(OUTMAN, "bad -k value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            UseKSVD=true;
            UseApproxKSVD =true;
            break;
            case 'm':
            if (sscanf(OptArg,"%s", nameInputMetric) != 1){
                fprintf(OUTMAN, "bad -m value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            UseInputMetric =true;
            break;
            case 'n':
            if(sscanf(OptArg,"%zu", &NAtoms) !=1) {
                fprintf(OUTMAN, "bad -n value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            break;
            case 'p':
            if (sscanf(OptArg,"%s", nameMeanPatch) != 1){
                fprintf(OUTMAN, "bad -p value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            CenterTraining =true;
            break;
            case 'r':
            if(sscanf(OptArg,"%lf", &InitDensity) !=1) {
                fprintf(OUTMAN, "bad -r value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            useSparseInit=true;
            break;
            case 's':
            if (sscanf(OptArg,"%d", &SparsityTarget) != 1){
                fprintf(OUTMAN, "bad -s value: %s\n", OptArg);
                exit(EXIT_FAILURE);
            }
            break;
            case 't':
            #ifdef _OPENMP
                if (sscanf(OptArg,"%d:%d",&outer_loop_threads,&inner_loop_threads) < 2){
                    printf( "Error: syntaxe is outer_loop_threads:inner_loop_threads ... \n");
                    exit(EXIT_FAILURE);
                }
                OMPThreadHandling::autoSetThreads(outer_loop_threads,inner_loop_threads,outer_loop_threads*inner_loop_threads);
            #endif
             break;

            case 'C': CenterPatches = true; break;
            case 'K': UseKSVD = true; break;
            case 'M':
                if (sscanf(OptArg,"%lf", &minimalCorrelation) != 1){
                    fprintf(OUTMAN, "bad -M value: %s\n", OptArg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'N':
                if (sscanf(OptArg,"%d", &minimalSparsity) != 1){
                    fprintf(OUTMAN, "bad -S value: %s\n", OptArg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'R': useRandomTrainInit = true; break;
            case 'T': timing = true; break;
            case 'V': verbose = true; break;
            case '?': usage(argv);
        }
    }
    if (OptInd < argc) strcpy(nameVectIn, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(nameDicoOut, argv[OptInd++]);
    else usage(argv);

    /* make sure there are not too many parameters */
    if (OptInd < argc){
        fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
        exit(EXIT_FAILURE);
    }
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    /* Get command line arguments, open input file(s) if necessary */
    transinit(argc, argv);
    if(verbose){
        std::cout << "Test Dico in = " << nameDicoIn << std::endl;
        std::cout << "Training Data in = " << nameVectIn << std::endl;
        std::cout << "Output Dico in = " << nameDicoOut << std::endl;
        if (UseSparsity) std::cout << "Sparsity target = " <<
                  minimalSparsity << "-" << SparsityTarget << std::endl;

        std::cout << "Use Error target " << std::endl;
    }
    if (verbose) std::cout << "Set the RNG"<< std::endl;

    //Get Output prefix
    int len_substr;
    char *pos_substr,prefix[1024],tempname[1024];
    pos_substr=strrchr(nameDicoOut,'.');
    if(pos_substr==NULL) pos_substr=strrchr(nameDicoOut,'\0');
    len_substr=pos_substr-nameDicoOut;
    strncpy(prefix, nameDicoOut,len_substr);
    prefix[len_substr]='\0';
    if(verbose) printf("prefix=%s\n",prefix);

    struct timeval start,end,diff;
    gettimeofday(&start, (struct timezone *) NULL);

    //Init Random Generator
    gslRNG* RNGenerator = new gslRNG();
    //RNGenerator->createSeed();
    RNGenerator->setSeed((size_t) GeneratorSeed);
    if (verbose) std::cout << "RNG seed="<< RNGenerator->getSeed()<<std::endl;

    //Get Training Data
    if (verbose) std::cout << "Reading the data"<< std::endl;
    dblarray *Dico, Data, TabPatch, SigBlock;
    fits_read_dblarr(nameVectIn, Data);
    size_t Nx = Data.nx();//Number of pixels
    size_t Ny = Data.ny();//Number of patches
    size_t Nz = Data.nz();
    if (verbose){
        if (Data.naxis() == 2)
        std::cout << "Training set: " << Ny << " samples of " << Nx
                                                        << " pixels"<<std::endl;
        else std::cout << "Training set: Nx = " << Nx << " Ny = " << Ny
                                                << " Nz = " << Nz << std::endl;
    }

    //Initialize Sparse Coder
    SparseCoding *SparseCoder = new OMP();

    //Set Input Structure
    std::vector<size_t> ChannelDim;
    ChannelDim.push_back(Nx);
    std::vector<Manifold> _manifold;
    _manifold.push_back(Manifold::Rn);//Euclidean
    DataStruc* DataStructure=new DataStruc(&Data,ChannelDim,_manifold,true);

    //Get Initial Dictionary
    DictionaryInit* InitDico;
    if(verbose) printf("useRandomTrainInit=%d\n", useRandomTrainInit?1:0);

    if(useRandomTrainInit){
        InitDico=new RandomDictInit(NAtoms,&Data);
        if (verbose)
        std::cout << "Set Initial dictionary from random training samples. "
                                                                    <<std::endl;
    } else {
        InitDensity=(InitDensity<Nx/(10.*Ny))? Nx/(10.*Ny): InitDensity;
        if (verbose)
        if(verbose) std::cout<<"Random sparse init dico with density: "
                                                    << InitDensity <<std::endl;
        InitDico=new SparseDictInit(InitDensity,NAtoms, Nx);
    }
    InitDico->setRNG(RNGenerator);

    if(ReadInitDico) {
        if(verbose)
        if(verbose) std::cout<<"Read Initial dico in: "<<nameDicoIn<<std::endl;
        Dico=InitDico->readDicoFits(nameDicoIn);
        sprintf(tempname,"%s_read_initdico.fits",prefix);
    } else {
        Dico=InitDico->computeInit();
        sprintf(tempname,"%s_%s_%s_initdico.fits",prefix,
        RNGenerator->getSignature().c_str(),InitDico->getSignature().c_str());
        fits_write_dblarr(tempname, *Dico);
    }
    delete InitDico;

    NAtoms = Dico->nx();//Number of Atoms
    if (verbose) {
        if (Dico->naxis() == 2)
        std::cout << "Initial dictionary: " << NAtoms << " atoms of " << Nx
                                                        << " pixels"<<std::endl;
        else std::cout << "Initial dictionary: Nx = " << Nx << " Ny = " << NAtoms
                                                << " Nz = " << Dico->nz() << endl;
    }
    if(NAtoms ==0) {
        std::cout << "NAtoms should be >0: " << NAtoms <<std::endl;
        std::cout << "Bad initialization. Use -n or different input?" <<std::endl;
        exit(EXIT_FAILURE);
    }

    if (SparsityTarget<0) SparsityTarget = NAtoms;

    SparseCoder->setErrorTarget(ErrorTarget);
    SparseCoder->setSparsityTarget(SparsityTarget);
    SparseCoder->setVerbose(verbose);
    //Initialize Matrix decomposition
    #if defined _ARMA_DEC
        MatrixDecomposition* ArmaDec = new ArmadilloDecomposition();
    #elif defined _ATLAS_DEC
        MatrixDecomposition* ArmaDec = new AtlasDecomposition();
    #else
        MatrixDecomposition* ArmaDec = new SapDecomposition();
    #endif
    SparseCoder->setMatDecomp(*ArmaDec);
    dblarray *InputMetric;
    if(UseInputMetric){
        SparseCoder->setDictionary(*Dico);
        InputMetric= SparseCoder->readMetric(nameInputMetric);
        SparseCoder->setInputMetric(*InputMetric);
    }
    //Compute the mean of input patches in case we want to add it back for
    // sparse approximation
    dblarray* meanArray;
    if(SaveApprox){
        if(CenterPatches) {
            if(verbose) std::cout<<"REMOVE MEAN PATCH"<<std::endl;
            meanArray = new dblarray();
            ArmaDec->getMeanAlongDim(Data, 0, meanArray);
        }
    }

    //Initialize Dictionary Learning
    SparseVectorDL *DictLearning;
    if(UseKSVD){
        DictLearning = new SparseVectorKSVD(*SparseCoder,UseApproxKSVD, KSVDit);
        if(verbose) std::cout<<"USE KSVD"<<std::endl;
    }  else DictLearning = new SparseVectorMOD(*SparseCoder,False);
    if(verbose) std::cout<<"SCODING="<<SparseCoder->getSignature()<<std::endl;
    if(DictLearning->getCoder() ==nullptr) std::cout<<"CODER IS NULL !"<<std::endl;
    DictLearning->setMatDecomp(*ArmaDec);
    DictLearning->setTrainingSetCentering(CenterTraining);
    DictLearning->setPatchCentering(CenterPatches);
    DictLearning->setRNG(RNGenerator);
    DictLearning->setDictionary(*Dico);
    DictLearning->setTraining(*DataStructure);
    DictLearning->setVerbose(verbose);
    if(NbIter>0) DictLearning->setNIter(NbIter);
    else DictLearning->setNIter(10);
    if (verbose) {
        std::cout << "DL " <<DictLearning->getSignature()<<": use "<<
                    DictLearning->getNAtoms() << " atoms of " <<
                    DictLearning->getNPix() << " pixels learned with "<<
                    DictLearning->getNTrain()<<" samples" <<std::endl;
    }

    //Finally, LEARN !
    DictLearning->learning();

    gettimeofday(&end, (struct timezone *) NULL);
    timersub(&end,&start,&diff);
    if(timing)
        std::cout << "Time for " <<DictLearning->getSignature() <<":"<<
                    (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;

    //Save Learned Dictionary
    fits_write_dblarr (nameDicoOut, *(DictLearning->getDictionary()));
    dblarray *code, *approx;
    if(SaveCode) {
        code=DictLearning->computeCode(true);
        fits_write_dblarr(nameCodeOut, *(code));
    }
    if(SaveApprox){
        approx=DictLearning->getApprox(!SaveCode);
        if(CenterPatches) {
            ArmaDec->addAlongDim(*approx, 0, *meanArray,1.0);
            delete meanArray;
        }
        fits_write_dblarr(nameApproxOut, *(approx));
    }

    if(verbose) printf("delete Data \n");
    delete DataStructure;
    Data.free();
    if(verbose) printf("delete DictLearning \n");
    delete DictLearning;
    if(verbose) printf("delete ArmaDec \n");
    delete ArmaDec;
    if(UseInputMetric){
        if(verbose) printf("delete input metric \n");
        delete InputMetric;
    }
    if(verbose) printf("delete SparseCoder \n");
    delete SparseCoder;
    if(verbose) printf("delete Dico \n");
    delete Dico;
    if(verbose) printf("delete RNGenerator \n");
    delete RNGenerator;
    MyExit::exit(EXIT_SUCCESS);
}

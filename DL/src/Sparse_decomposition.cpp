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
 **    File:  Sparse_decomposition.cpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION Main program to perform sparse decomposition given a dictionary
 **    -----------
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
#include "SparseCoding.h"
#include "OMPCoding.h"
#ifdef _OPENMP
#include "OpenMPInterface.h"
#endif

//#define _DEBUG_COST
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

//Code for input data, with Npatches rows and Npix columns, stored in row
//major order.
// Beware of dblarray conventions: Index katom,ktrain of CodeSet accessed as:
// Codeset(katom,kpatch)=buffer[katom+ kpatch*Natoms];
// This allows to access fastly to all coefs of a given patch (e.g. for sparse
// coding)
char nameCodeOut[1024];

//optional output approximation for patches, same convention as input patches.
char nameApproxOut[1024];
//optional mean patch from training set.
char nameMeanInput[1024];
//optional input metric to use for sparse coding atom selection.
char nameInputMetric[1024];
//optional input weights to use for sparse coding selection and approximation
char nameInputWeights[1024];
//option reference for S1 log-map computation
char nameInputReferenceS1[1024];

bool verbose=false;
bool timing=false;
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, const char *opts);
double EpsCompare=1e-15;//Expected accuracy
size_t bufferSize=0;

double ErrorTarget = -1.;
size_t NAtoms=0;
bool UseSparsity=true;
int SparsityTarget=-1;
bool SaveApprox=false;

bool CenterPatches=false;
bool CenterTraining=false;
bool UseInputMetric =false;
bool UseInputWeights=false;
bool SaveCode=True;

/*********************************************************************/
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options nameDicoIn.fits", argv[0]);
    fprintf(OUTMAN, " nameVectIn.fits nameCodeOut.fits \n\n");
    fprintf(OUTMAN, "   where options =  \n");
    fprintf(OUTMAN, "      [-e ErrorTarget, TargetSparsity ignored]\n");
    fprintf(OUTMAN, "      [-s TargetSparsity]\n");
    fprintf(OUTMAN, "      [-m name of fits file containing input metric for selection of atoms]\n");
    fprintf(OUTMAN, "      [-w name of fits file containing input metric for selection and approximation]\n");
    fprintf(OUTMAN, "      [-p Name of optional mean input to subtract\n");
    fprintf(OUTMAN, "          from each input]\n");
    fprintf(OUTMAN, "      [-a Name for training set approx to save]\n");
    fprintf(OUTMAN, "      [-C center training patches (mean of each patch=0)]\n");
    fprintf(OUTMAN, "      [-N do not save the final code]\n");
    fprintf(OUTMAN, "      [-T Timing on]\n");
    fprintf(OUTMAN, "      [-V verbose]\n");
    MyExit::exit(EXIT_SUCCESS);
}

/*********************************************************************/
/* GET COMMAND LINE ARGUMENTS */
static void transinit(int argc, char *argv[])
{
    int c,outer_loop_threads,inner_loop_threads;
    /* get options */
    while ((c = GetOpt(argc,argv,"a:b:e:m:p:s:t:z:w:CHTVXZN")) != -1){
        switch (c){
            case 'a':
            if (sscanf(OptArg,"%s", nameApproxOut) != 1){
                fprintf(OUTMAN, "bad -a value: %s\n", OptArg);
                MyExit::exit(EXIT_FAILURE);
            }
            SaveApprox =true;
            break;
            case 'b':
            if (sscanf(OptArg,"%lu", &bufferSize) != 1){
                fprintf(OUTMAN, "bad -b value: %s\n", OptArg);
                MyExit::exit(EXIT_FAILURE);
            }
            break;
            case 'e':
            if(sscanf(OptArg,"%lf", &ErrorTarget) !=1){
                fprintf(OUTMAN, "bad -e value: %s\n", OptArg);
                MyExit::exit(EXIT_FAILURE);
            }
            break;
            case 'm':
            if (sscanf(OptArg,"%s", nameInputMetric) != 1){
                fprintf(OUTMAN, "bad -m value: %s\n", OptArg);
                MyExit::exit(EXIT_FAILURE);
            }
            UseInputMetric =true;
            break;
            case 't':
            #ifdef _OPENMP
                std::cout<<"OPTARG:"<<OptArg<<std::endl;
                if(sscanf(OptArg,"%d:%d",&outer_loop_threads,\
                                                    &inner_loop_threads) < 2){
                    printf("Error in setting threads \n");
                    exit(EXIT_FAILURE);
                }
                OMPThreadHandling::autoSetThreads(outer_loop_threads,
                    inner_loop_threads,outer_loop_threads*inner_loop_threads);
            #endif
             break;
            case 'p':
            if (sscanf(OptArg,"%s", nameMeanInput) != 1){
                fprintf(OUTMAN, "bad -p value: %s\n", OptArg);
                MyExit::exit(EXIT_FAILURE);
            }
            CenterTraining =true;
            break;
            case 's':
            if (sscanf(OptArg,"%d", &SparsityTarget) != 1){
                fprintf(OUTMAN, "bad -s value: %s\n", OptArg);
                MyExit::exit(EXIT_FAILURE);
            }
            break;
            case 'w':
            if (sscanf(OptArg,"%s", nameInputWeights) != 1){
                fprintf(OUTMAN, "bad -w value: %s\n", OptArg);
                MyExit::exit(EXIT_FAILURE);
            }
            UseInputWeights =true;
            break;
            case 'C': CenterPatches = true; break;
            case 'N': SaveCode = false; break;
            case 'T': timing = true; break;
            break;
            case 'V': verbose = true; break;
            case '?': usage(argv);
        }
    }
    if (OptInd < argc) strcpy(nameDicoIn, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(nameVectIn, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(nameCodeOut, argv[OptInd++]);
    else usage(argv);

    /* make sure there are not too many parameters */
    if (OptInd < argc){
        fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
        MyExit::exit(EXIT_FAILURE);
    }
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    /* Get command line arguments, open input file(s) if necessary */
    transinit(argc, argv);
    if(verbose){
        std::cout << "Input Dico in = " << nameDicoIn << std::endl;
        std::cout << "Data 2 sparsify in = " << nameVectIn << std::endl;
        std::cout << "Output Code in = " << nameCodeOut << std::endl;
        if (UseSparsity) std::cout << "Use Sparsity target = " << SparsityTarget
                                                                    << std::endl;
        else std::cout << "Use Error target " << std::endl;
    }
    if (verbose) std::cout << "Set the RNG"<< std::endl;

    //Get Output prefix
    int len_substr;
    char *pos_substr,prefix[1024];
    pos_substr=strrchr(nameCodeOut,'.');
    if(pos_substr==NULL) pos_substr=strrchr(nameCodeOut,'\0');
    len_substr=pos_substr-nameCodeOut;
    strncpy(prefix, nameCodeOut,len_substr);
    prefix[len_substr]='\0';
    if(verbose) printf("prefix=%s\n",prefix);

    struct timeval start,end,diff;
    gettimeofday(&start, (struct timezone *) NULL);
    OMPThreadHandling::getStatus();

    //Get Input Data
    if (verbose) std::cout << "Reading the data"<< std::endl;
    dblarray Dico, Data, TabPatch, SigBlock;


    fits_read_dblarr(nameVectIn, Data);//Should be vectorized at one time
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
    //Get Dictionary
    if(verbose) printf("DATA READ DICO %s\n",nameDicoIn);
    fits_read_dblarr (nameDicoIn, Dico);
    NAtoms = Dico.nx();//Number of Atoms
    Nx = Dico.ny();//Vectorize Patch
    if (verbose) {
        if (Dico.naxis() == 2)
            std::cout << "Input dictionary: " << NAtoms << " atoms of " << Nx
                                                        << " pixels"<<std::endl;
        else std::cout << "Input dictionary: Nx = " << Nx << " Ny = " << NAtoms
                                                << " Nz = " << Dico.nz() << endl;
    }
    if(NAtoms <=0) {
        std::cout << "NAtoms should be >0: " << NAtoms <<std::endl;
        std::cout << "Bad Dictionary. Correct path ?" <<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    if(!SaveCode && !SaveApprox){
        std::cout << "No code nor approximation to save " <<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }

    //Initialize Sparse Coder
        //Initialize Sparse Coder
    SparseCoding *SparseCoder = new OMP();
    if (SparsityTarget<0) SparsityTarget = NAtoms;
    SparseCoder->setDictionary(Dico);
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
        InputMetric= SparseCoder->readMetric(nameInputMetric);
        SparseCoder->setInputMetric(*InputMetric);
    }
    SparseCoder->setVerbose(verbose);
    dblarray* meanArray;
    if(CenterPatches) {
        if(verbose) std::cout<<"REMOVE MEAN PATCH"<<std::endl;
        meanArray = new dblarray();
        ArmaDec->centerRow2DData(Data,meanArray);
    } else if(CenterTraining){
        if(verbose) std::cout<<"Center Training"<<std::endl;
        meanArray = new dblarray();
        fits_read_dblarr(nameMeanInput, *meanArray);
        ArmaDec->subCol2DData(Data, *meanArray);
        delete meanArray;
    }

    if (verbose) {
        std::cout << "SparseCoding " <<SparseCoder->getSignature()<<": use "<<
                        SparseCoder->getNbAtoms() << " atoms of " <<
                        SparseCoder->getNbpixels() << " pixels to code "<<
                        Data.ny()<<" samples" <<std::endl;
    }

    //Set Input Structure
    std::vector<size_t> ChannelDim;
    ChannelDim.push_back(Nx);
    std::vector<Manifold> _manifold;
    _manifold.push_back(Manifold::Rn);//Euclidean

    DataStruc* DataStructure;
    dblarray *WeightedInput;
    dblarray* SparseCode;
    if(UseInputWeights){
            dblarray *InputWeights= SparseCoder->readMetric(nameInputWeights);
            dblarray* TempData=&Data;
            WeightedInput=SparseCoder->applySqrtWeights(TempData,*InputWeights);
            delete InputWeights;
            Data.free();
            DataStructure= new DataStruc(WeightedInput,ChannelDim,_manifold,true);
    } else DataStructure= new DataStruc(&Data,ChannelDim,_manifold,true);

    //Finally, get Sparse Code !
    if((SaveCode)||(bufferSize==0)) {//Need to bufferize this also
        if(verbose) std::cout<<"SAVE CODE"<<std::endl;
        SparseCode=SparseCoder->computeCode(DataStructure);
        gettimeofday(&end, (struct timezone *) NULL);
        timersub(&end,&start,&diff);
        if(timing)
            std::cout << "Time to get code: " <<SparseCoder->getSignature() <<":"<<
                    (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
        //Save Sparse Code Dictionary
        fits_write_dblarr(nameCodeOut, *(SparseCode));
	    if((!SaveApprox)||((SaveApprox)&&(bufferSize>0))){
	        printf("delete SparseCode \n");
            delete SparseCode;
        }
	}

    if(SaveApprox){
        dblarray* Approx;
	    DataStruc* NewStruc;
        if(bufferSize>0){//Need to bufferize reading/writing ideally
            size_t Npatches=DataStructure->getNPatches();
            size_t TwoDPatchSize=DataStructure->get2DPatchSize();
            size_t offsetPatches=0;
	        Approx=new dblarray();
            Approx->alloc(TwoDPatchSize,Npatches);
            dblarray* ApproxTemp=new dblarray();
            double* OffsetStart=Approx->buffer();
            while(offsetPatches+bufferSize<Npatches){
    		    if(verbose) std::cout<<"OFFSET="<<offsetPatches<<std::endl;
                NewStruc=DataStructure->extractContiguousPatches(offsetPatches,
                                                                    bufferSize);
                SparseCode=SparseCoder->computeCode(NewStruc);
                ApproxTemp->alloc(OffsetStart,TwoDPatchSize,bufferSize,0);
                SparseCoder->computeApprox(SparseCode,ApproxTemp,NewStruc);
                OffsetStart+=bufferSize*TwoDPatchSize;
                offsetPatches+=bufferSize;
                delete NewStruc;
                delete SparseCode;
            }
            size_t remPatches=Npatches-offsetPatches;
	        if(remPatches>0) {
            	NewStruc=DataStructure->extractContiguousPatches(offsetPatches,
                                                                remPatches);
            	SparseCode=SparseCoder->computeCode(NewStruc);
            	ApproxTemp->alloc(OffsetStart,TwoDPatchSize,remPatches,0);
            	SparseCoder->computeApprox(SparseCode,ApproxTemp,NewStruc);
            	delete SparseCode;
	        }
            delete NewStruc;
            delete ApproxTemp;
            if(verbose) std::cout<<"ADD BACK MEAN"<<std::endl;
            if(CenterPatches) {
                ArmaDec->addAlongDim(*Approx, 0, *meanArray,1.0);
                delete meanArray;
            }
            fits_write_dblarr(nameApproxOut, *(Approx));
            delete Approx;
        } else {
            dblarray* approx;
            approx=SparseCoder->computeApprox(SparseCode,DataStructure);
            if(verbose) std::cout<<"ADD BACK MEAN"<<std::endl;
            if(CenterPatches) {
                ArmaDec->addAlongDim(*approx, 0, *meanArray,1.0);
                delete meanArray;
            }
            fits_write_dblarr(nameApproxOut, *(approx));
            delete approx;
    	}
		delete SparseCode;
    }

    delete DataStructure;
    if(UseInputWeights){
        if(verbose) printf("delete weighted input \n");
        delete WeightedInput;
    } else Data.free();
    printf("delete ArmaDec \n");
    delete ArmaDec;
    if(UseInputMetric){
        if(verbose) printf("delete input metric \n");
        delete InputMetric;
    }
    if(verbose) printf("delete OMPCoding \n");
    delete SparseCoder;
    MyExit::exit(EXIT_SUCCESS);
}

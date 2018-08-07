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
 **    File:  Dl_test_matrix_decomposition.cpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION  Unit tests to check matrix decompositions
 **    -----------
 **
 **    Usage: DL_test_matrix_decomposition options output
 **
 ******************************************************************************/
#include "MatrixDecomposition.h"
#include <cmath>
#include "Array.h"
#include "NR.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM1D_Block.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iomanip>
#include <sys/time.h>

/****************************************************************************/
char nameTest[1024];  // 2D array, containing fits
char nameOutputTest[1024];  // 2D array
bool testInput=false;
bool verbose=false;
bool testChol=false;
bool testSVD=false;
bool testInvert=false;
bool testPInvert=false;
bool timing=false;
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char *const*argv, const char *opts);
double EpsCompare=1e-15;//Expected accuracy

/*********************************************************************/
static void usage(char *argv[])
{
   fprintf(OUTMAN, "Usage: %s options nameOutputTest.fits \n\n", argv[0]);
   fprintf(OUTMAN, "   where options =  \n");
   fprintf(OUTMAN, "        [-I flag to test inversion]\n");
   fprintf(OUTMAN, "        [-C flag to test cholesky decomposition]\n");
   fprintf(OUTMAN, "        [-P flag to test pseudo-inverse]\n");
   fprintf(OUTMAN, "        [-S flag to test SVD]\n");
   fprintf(OUTMAN, "        [-T Timing on]\n");
   fprintf(OUTMAN, "        [-t nameTest]\n");
   fprintf(OUTMAN, "        [-V verbose]\n");
   exit(EXIT_SUCCESS);
}
/*********************************************************************/
/* GET COMMAND LINE ARGUMENTS */
static void transinit(int argc, char *argv[])
{
   int c;
   /* get options */
   while ((c = GetOpt(argc,argv,"CIPSTt:V")) != -1){
      switch (c){
         case 'C': testChol = true; break;
         case 'I': testInvert = true; break;
         case 'P': testPInvert = true; break;
         case 'S': testSVD = true; break;
         case 'T': timing = true; break;
         break;
         case 't':
            if (sscanf(OptArg,"%s", nameTest) != 1){
               fprintf(OUTMAN, "bad -s value: %s\n", OptArg);
               exit(EXIT_FAILURE);
            }
            testInput =true;
         break;
         case 'V': verbose = true; break;
         case '?': usage(argv);
      }
   }
   if (OptInd < argc) strcpy(nameOutputTest, argv[OptInd++]);
   else {
      if(testInput) usage(argv);
      else sprintf(nameOutputTest,"default_output.fits");

   }
   /* make sure there are not too many parameters */
   if (OptInd < argc){
      fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
      exit(EXIT_FAILURE);
   }
}


/****************************************************************************/
dblarray* init_simple_square_array() {
   double x[] = {1,4,2,1,0,2,2,2,0};
   dblarray* testSquareArray= new dblarray(3,3,0);
   double *buf= testSquareArray->buffer();
   for(int k=0;k<9;++k,++buf) *buf=x[k];
   return testSquareArray;
}

bool check_inverse_simple_square_array(dblarray &testArray) {
   double xinv[]={-0.25,0.25,0.5,0.25,-0.25,0.,0.125,0.375,-0.25};
   double *buf= testArray.buffer();
   bool check=true;
   for(int k=0;k<9;++k,++buf){
      if(fabs((*buf)-xinv[k])> EpsCompare){
         check=false;
         std::cout<<"Expecting inverse " <<std::setprecision(15)<<std::fixed<<
                                       xinv[k] << " , got " << *buf<<std::endl;
      }
   }
   return check;
   //inv:[[-0.5,0.5,0.25],[0.5,-0.5,0.25],[0.25,0.25,-0.125]]
}

dblarray* init_simple_rect_array() {
   double x[] = {1,1,-sqrt(3)/4.,sqrt(3)/4.,0.25,-0.25};
   dblarray* testRectArray= new dblarray(2,3,0);//2 columns, 3 lines
   double *buf= testRectArray->buffer();
   for(int k=0;k<6;++k,++buf) *buf=x[k];
   return testRectArray;
}

dblarray* init_square_array_svd() {
   double x[] = {42/16.0,-22/16.0,-sqrt(6)/8.0,-22/16.0,42/16.0,-sqrt(6)/8.0,
                                            -sqrt(6)/8.0,-sqrt(6)/8.0,7/4.0};
   dblarray* testSqrArray= new dblarray(3,3,0);//3 columns, 3 lines
   double *buf= testSqrArray->buffer();
   for(int k=0;k<9;++k,++buf) *buf=x[k];
   return testSqrArray;
}

bool check_svd_simple_rect_array(dblarray &U, dblarray &S, dblarray &Vt) {
   double svdU[]={1.,0.,0.,sqrt(3.)/2.,0,-0.5};
   double svdS[]={sqrt(2.),sqrt(2.)/2.};
   double svdVt[]={sqrt(2)/2,sqrt(2)/2,-sqrt(2)/2,sqrt(2)/2};
   double *bufVt=Vt.buffer(),*bufS=S.buffer(), *bufU=U.buffer();
   bool check=true;

   if(verbose){
      std::cout << "Usize:"<<U.nx()<<","<<U.ny()<<std::endl;
      std::cout << "Ssize:"<<S.nx()<<","<<S.ny()<<std::endl;
      std::cout << "Vtsize:"<<Vt.nx()<<","<<Vt.ny()<<std::endl;
   }
   //Lift some arbitrary choices
   if(bufS[0]<0) {
      bufS[0]*=-1;
      bufVt[0]*=-1;
      bufVt[1]*=-1;
   }
   if(bufS[1]<0) {
      bufVt[2]*=-1;
      bufVt[3]*=-1;
   }
   if(bufVt[0] < 0) {
      bufVt[0]*=-1;
      bufVt[1]*=-1;
      bufU[0]*=-1;
   }
   if(bufVt[2] > 0) {
      bufVt[2]*=-1;
      bufVt[3]*=-1;
      bufU[3]*=-1;
      bufU[5]*=-1;
   }
   //First check Vt
   for(int k=0;k<4;++k,++bufVt){
      if(fabs((*bufVt)-svdVt[k])> EpsCompare){
         check=false;
         std::cout<<"Expecting Vt[" <<std::setprecision(15)<<std::fixed<< k <<
                              "] "<< svdVt[k] << " , got " << *bufVt<<std::endl;
      }
   }
   //Check S
   for(int k=0;k<2;++k,++bufS){
      if(fabs((*bufS)-svdS[k])> EpsCompare){
         check=false;
         std::cout<<"Expecting S[" <<std::setprecision(15)<<std::fixed<< k <<
                              "] "<< svdS[k] << " , got " << *bufS <<std::endl;
      }
   }
   //Check U
   for(int k=0;k<6;++k,++bufU){
      if(fabs((*bufU)-svdU[k])> EpsCompare){
         check=false;
         std::cout<<"Expecting U["<<std::setprecision(15)<<std::fixed<< k <<
                                "] "<< svdU[k] << " , got " << *bufU<<std::endl;
      }
   }
   return check;
}

bool check_pinverse_simple_rect_array(dblarray &testArray) {
   double xpinv[]={0.5,-sqrt(3.)/2.,0.5,0.5,sqrt(3.)/2.,-0.5};
   double *buf= testArray.buffer();
   bool check=true;
   for(int k=0;k<6;++k,++buf){
      if(fabs((*buf)-xpinv[k])> EpsCompare){
         check=false;
         std::cout<<"Expecting Pseudo-inverse " <<std::setprecision(15)<<std::fixed<<
                                       xpinv[k] << " , got " << *buf<<std::endl;
      }
   }
   return check;
}

bool check_sqrtSVD_simple_array(dblarray &testArray) {
   double xsqrt[]={(11.0+sqrt(2.0))/8.0,(sqrt(2.0)-5.0)/8.0,(sqrt(6.0)-2*sqrt(3))/8.0,
        (sqrt(2.0)-5.0)/8.0,(11.0+sqrt(2.0))/8.0,(sqrt(6.0)-2*sqrt(3))/8.0,
        (sqrt(6.0)-2*sqrt(3))/8.0,(sqrt(6.0)-2*sqrt(3))/8.0,(1.0+3*sqrt(2))/4.0};

   double *buf= testArray.buffer();
   bool check=true;
   for(int k=0;k<9;++k,++buf){
      if(fabs((*buf)-xsqrt[k])> EpsCompare){
         check=false;
         std::cout<<"Expecting sqrt " <<std::setprecision(15)<<std::fixed<<
                                       xsqrt[k] << " , got " << *buf<<std::endl;
      }
   }
   return check;
}

dblarray* generate_random_array(size_t Nx,size_t Ny,size_t Nz,
                                                   double mean, double sigma){
   gsl_rng *rng;
   dblarray *outRand=new dblarray(Nx,Ny,Nz);
   const gsl_rng_type * T;
   T = gsl_rng_default;
   rng = gsl_rng_alloc (T);

   double *buffer=outRand->buffer();
   for(size_t kel=0;kel<(size_t)outRand->n_elem();++kel,++buffer)
      *buffer=(gsl_ran_gaussian(rng,sigma)+mean);
   gsl_rng_free(rng);
   return outRand;
}

void check_scale_centering_routines(MatrixDecomposition* MatDec){
   dblarray *randColArray=generate_random_array(1000,5000,0,2.,4.);
   dblarray *randRowArray=generate_random_array(5000,1000,0,2.,5.);
   dblarray copyColArray= *randColArray;
   dblarray copyColArray2= *randColArray;
   dblarray copyRowArray= *randRowArray;
   dblarray copyRowArray2= *randRowArray;

   dblarray subArray;
   double mean;
   struct timeval start,end,diff;
   //Test Centering routine for Column Data
   gettimeofday(&start, (struct timezone *) NULL);
   for(size_t kx=0;kx<(size_t)randColArray->nx();++kx) {//1000 vecof 5000 el
      mean=0.;
      for(size_t ky=0;ky<(size_t)randColArray->ny();++ky)
                                                     mean+= copyColArray(kx,ky);
      mean/=randColArray->ny();
      for(size_t ky=0;ky<(size_t)randColArray->ny();++ky)
                                                      copyColArray(kx,ky)-=mean;
   }
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   std::cout << "\t Time to manually center Col data :"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   gettimeofday(&start, (struct timezone *) NULL);
   MatDec->centerCol2DData(*randColArray);
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   std::cout << "\t Time to center Col data for " <<MatDec->getSignature() <<":"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   subArray= copyColArray-*randColArray;
   std::cout << "Max absolute difference for " <<MatDec->getSignature() <<":"<<
                 subArray.maxfabs()<<std::endl;
   //Test Centering routine for Row Data
   gettimeofday(&start, (struct timezone *) NULL);
   for(size_t ky=0;ky<(size_t)randRowArray->ny();++ky) {//1000 vec of 5000 els
      mean=0.;
      for(size_t kx=0;kx<(size_t)randRowArray->nx();++kx)
                                                    mean+= copyRowArray(kx,ky);
      mean/=randRowArray->nx();
      for(size_t kx=0;kx<(size_t)randRowArray->nx();++kx)
                                                      copyRowArray(kx,ky)-=mean;
   }
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   std::cout << "\t Time to manually center Row data :"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   gettimeofday(&start, (struct timezone *) NULL);
   MatDec->centerRow2DData(*randRowArray);
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   std::cout << "\t Time to center Row data for " <<MatDec->getSignature() <<":"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   subArray= copyRowArray-* randRowArray;
   std::cout << "Max absolute difference for " <<MatDec->getSignature() <<":"<<
                 subArray.maxfabs()<<std::endl;

   //Test scaling routine for Column Data
   gettimeofday(&start, (struct timezone *) NULL);
   for(size_t kx=0;kx<(size_t)randColArray->nx();++kx) {//1000 vec of 5000 els
      mean=0.;
      for(size_t ky=0;ky<(size_t)randColArray->ny();++ky)
                                  mean+=copyColArray(kx,ky)*copyColArray(kx,ky);
      for(size_t ky=0;ky<(size_t)randColArray->ny();++ky)
                                     if(mean>1) copyColArray(kx,ky)/=sqrt(mean);
   }
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   std::cout << "\t Time to manually scale Col data :"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   gettimeofday(&start, (struct timezone *) NULL);
   MatDec->l2ProjCol2DData(*randColArray);
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   std::cout << "\t Time to scale data Col for " <<MatDec->getSignature() <<":"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   subArray.init(-1);
   subArray= copyColArray-*randColArray;
   std::cout << "Max absolute difference for " <<MatDec->getSignature() <<":"<<
                 subArray.maxfabs()<<std::endl;
   //Test scaling routine for Row Data
   gettimeofday(&start, (struct timezone *) NULL);
   for(size_t ky=0;ky<(size_t)randRowArray->ny();++ky) {//1000 vec of 5000 els
      mean=0.;
      for(size_t kx=0;kx<(size_t)randRowArray->nx();++kx)
                                  mean+=copyRowArray(kx,ky)* copyRowArray(kx,ky);
      for(size_t kx=0;kx<(size_t)randRowArray->nx();++kx)
                                     if(mean>1) copyRowArray(kx,ky)/=sqrt(mean);
   }
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   std::cout << "\t Time to manually scale Row data :"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   gettimeofday(&start, (struct timezone *) NULL);
   MatDec->l2ProjRow2DData(*randRowArray);
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   std::cout << "\t Time to scale Row data for " <<MatDec->getSignature() <<":"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   subArray.init(-1);
   subArray= copyRowArray-*randRowArray;
   std::cout << "Max absolute difference for " <<MatDec->getSignature() <<":"<<
                 subArray.maxfabs()<<std::endl;

   //Test Normalization routine for Column Data
   gettimeofday(&start, (struct timezone *) NULL);
   MatDec->l2NormalizeCol2DData(copyColArray2);
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   std::cout << "\t Time to norm col data for " <<MatDec->getSignature() <<":"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   subArray.init(-1);
   subArray= copyColArray2-copyColArray;
   std::cout << "Max absolute difference for " <<MatDec->getSignature() <<":"<<
                 subArray.maxfabs()<<std::endl;
   //Test Normalization routine for Row Data
   gettimeofday(&start, (struct timezone *) NULL);
   MatDec->l2NormalizeRow2DData(copyRowArray2);
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   std::cout << "\t Time to norm row data for " <<MatDec->getSignature() <<":"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   subArray.init(-1);
   subArray= copyRowArray2-copyRowArray;
   std::cout << "Max absolute difference for " <<MatDec->getSignature() <<":"<<
                 subArray.maxfabs()<<std::endl;
   delete randColArray;
   delete randRowArray;

}


void check_invert(MatrixDecomposition* MatDec,dblarray &Data,char* prefix,
                   const char* namespec,bool timer=false){
   struct timeval start,end,diff;
   char tempname[1024];
   gettimeofday(&start, (struct timezone *) NULL);
   dblarray *MatrixInvert= MatDec->invert(Data);
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   if(timer)
      std::cout << "\t Time " <<MatDec->getSignature() <<":"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;

   sprintf(tempname,"%s_%s_%s.fits",prefix, namespec,
                                             MatDec->getSignature().c_str());
   fits_write_dblarr(tempname,*MatrixInvert);
   delete MatrixInvert;
}

void check_sqrtSVD(MatrixDecomposition* MatDec,dblarray &Data,char* prefix,
                   const char* namespec,bool timer=false){
   struct timeval start,end,diff;
   char tempname[1024];
   dblarray U,S,Vt;
   gettimeofday(&start, (struct timezone *) NULL);
   dblarray* sqrtMat=MatDec->sqrtSVD(Data,0.);
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   if(timer)
      std::cout << "\t Time " <<MatDec->getSignature() <<":"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;

   sprintf(tempname,"%s_%s_sqrtMat_%s.fits",prefix, namespec,
                                             MatDec->getSignature().c_str());
   fits_write_dblarr(tempname,*sqrtMat);
   dblarray* sqrtMat2=MatDec->matMult(*sqrtMat,*sqrtMat);
   delete sqrtMat;
   sprintf(tempname,"%s_%s_sqrtMat2_%s.fits",prefix, namespec,
                                             MatDec->getSignature().c_str());
   fits_write_dblarr(tempname,*sqrtMat2);
   delete sqrtMat2;
}



void check_svd(MatrixDecomposition* MatDec,dblarray &Data,char* prefix,
                   const char* namespec,bool timer=false){
   struct timeval start,end,diff;
   char tempname[1024];
   dblarray U,S,Vt;
   gettimeofday(&start, (struct timezone *) NULL);
   MatDec->svd(Data,U,S,Vt);
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   if(timer)
      std::cout << "Time " <<MatDec->getSignature() <<":"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;

   sprintf(tempname,"%s_%s_U_%s.fits",prefix, namespec,
                                             MatDec->getSignature().c_str());
   fits_write_dblarr(tempname,U);
   sprintf(tempname,"%s_%s_S_%s.fits",prefix, namespec,
                                             MatDec->getSignature().c_str());
   fits_write_dblarr(tempname,S);
   sprintf(tempname,"%s_%s_Vt_%s.fits",prefix, namespec,
                                             MatDec->getSignature().c_str());
   fits_write_dblarr(tempname,Vt);
}


void check_pinvert(MatrixDecomposition* MatDec,dblarray &Data,char* prefix,
                   const char* namespec,bool timer=false){
   struct timeval start,end,diff;
   char tempname[1024];
   gettimeofday(&start, (struct timezone *) NULL);
   dblarray *MatrixPInvert= MatDec->pinvert(Data);
   gettimeofday(&end, (struct timezone *) NULL);
   timersub(&end,&start,&diff);
   if(timer)
      std::cout << "Time " <<MatDec->getSignature() <<":"<<
                  (double) diff.tv_sec+(double)diff.tv_usec/1000000.<<std::endl;
   sprintf(tempname,"%s_%s_%s.fits",prefix, namespec,
                                                MatDec->getSignature().c_str());
   fits_write_dblarr(tempname,*MatrixPInvert);
   delete MatrixPInvert;
}



/*********************************************************************/

int main(int argc, char *argv[])
{
   /* Get command line arguments, open input file(s) if necessary */
   transinit(argc, argv);
   if(verbose){
      if(testInput) std::cout << "Test Matrix in = " << nameTest<< std::endl;
      std::cout << "Output Test in = " << nameOutputTest << std::endl;
   }

   //get root name for output files
   int len_substr;
   char *pos_substr,prefix[1024],tempname[1024];
   pos_substr=strrchr(nameOutputTest,'.');
   if(pos_substr==NULL) pos_substr=strrchr(nameOutputTest,'\0');
   len_substr=pos_substr-nameOutputTest;
   strncpy(prefix, nameOutputTest,len_substr);
   prefix[len_substr]='\0';
   if(verbose) printf("prefix=%s\n",prefix);

   MatrixDecomposition* SapDec=new SapDecomposition();
   #ifdef _ATLAS_DEC
      MatrixDecomposition* AtlasDec=new AtlasDecomposition();
   #endif
   #ifdef _ARMA_DEC
      MatrixDecomposition* ArmaDec=new ArmadilloDecomposition();
      MatrixDecomposition* ArmaDecStd=new ArmadilloDecomposition(STD);
   #endif
   dblarray *MatrixInvert,Data, *MatrixPInvert;

   check_scale_centering_routines(SapDec);
   #ifdef _ATLAS_DEC
      check_scale_centering_routines(AtlasDec);
   #endif
   #ifdef _ARMA_DEC
      check_scale_centering_routines(ArmaDec);
   #endif
   if(testInput){
      std::cout << "Read input matrix."<<std::endl;
      fits_read_dblarr(nameTest, Data);
   }
   int Nx = Data.nx(), Ny = Data.ny();
   if(verbose)
         std::cout << "Input matrix: [" << Nx <<","<<Ny<<"]"<<std::endl;

   if(testInvert){
      std::cout << "Test inversion "<<std::endl;
      dblarray *testSquareArray= init_simple_square_array();
      std::cout << "Test Decomposition "<<SapDec->getSignature()<<std::endl;
      MatrixInvert=SapDec->invert(*testSquareArray);
      check_inverse_simple_square_array(*MatrixInvert);
      delete MatrixInvert;
      #ifdef _ATLAS_DEC
         std::cout<< "Test Decomposition "<<AtlasDec->getSignature()<<std::endl;
         MatrixInvert= AtlasDec->invert(*testSquareArray);
         check_inverse_simple_square_array(*MatrixInvert);
         delete MatrixInvert;
      #endif
      #ifdef _ARMA_DEC
         std::cout << "Test Decomposition "<<ArmaDec->getSignature()<<std::endl;
         MatrixInvert= ArmaDec->invert(*testSquareArray);
         check_inverse_simple_square_array(*MatrixInvert);
         delete MatrixInvert;
      #endif
      delete testSquareArray;
      if(testInput){
         if(Nx!=Ny)
            std::cout << "Cannot perform inversion of non square matrix: [" <<
                                                  Nx <<","<<Ny<<"]"<<std::endl;
         else{
            std::string suffix=std::string("invert");
            #ifdef _ARMA_DEC
               check_invert(ArmaDec,Data,prefix, suffix.c_str(),timing);
            #endif
            #ifdef _ATLAS_DEC
               check_invert(AtlasDec,Data,prefix, suffix.c_str(),timing);
            #endif
            check_invert(SapDec,Data,prefix, suffix.c_str(),timing);
         }
      } else {
          std::string suffix=std::string("invert");
          dblarray *testSquareArray=generate_random_array(1000,1000,0,0,1);
          sprintf(tempname,"%s_input_inversion.fits",prefix);
          fits_write_dblarr(tempname,*testSquareArray);
          std::cout << "Test time to invert matrix: [" <<testSquareArray->nx()
                                  <<","<<testSquareArray->ny()<<"]"<<std::endl;
          std::cout << "Input: "<<tempname <<std::endl;
          #ifdef _ARMA_DEC
          check_invert(ArmaDec,*testSquareArray,prefix, suffix.c_str(),timing);
          #endif
          #ifdef _ATLAS_DEC
          check_invert(AtlasDec,*testSquareArray,prefix, suffix.c_str(),timing);
          #endif
          check_invert(SapDec,*testSquareArray,prefix, suffix.c_str(),timing);
         delete testSquareArray;
      }
   }

   if(testSVD) {
      dblarray U,S,Vt;
      std::cout << "Test SVD "<<std::endl;
      dblarray *testRectArray= init_simple_rect_array();
      std::cout << "Test Decomposition "<<SapDec->getSignature()<<std::endl;
      SapDec->svd(*testRectArray,U,S,Vt);
      check_svd_simple_rect_array(U,S,Vt);
      U.free();S.free();Vt.free();
      #ifdef _ATLAS_DEC
        std::cout << "Test Decomposition "<<AtlasDec->getSignature()<<std::endl;
         AtlasDec->svd(*testRectArray,U,S,Vt);
         check_svd_simple_rect_array(U,S,Vt);
         U.free();S.free();Vt.free();
      #endif
      #ifdef _ARMA_DEC
         std::cout << "Test Decomposition "<<ArmaDec->getSignature()<<std::endl;
         ArmaDec->svd(*testRectArray,U,S,Vt);
         check_svd_simple_rect_array(U,S,Vt);
         U.free();S.free();Vt.free();
         ArmaDecStd->svd(*testRectArray,U,S,Vt);
         check_svd_simple_rect_array(U,S,Vt);
         U.free();S.free();Vt.free();
      #endif
      delete testRectArray;
      if(testInput){
         std::string suffix=std::string("svd");
            #ifdef _ARMA_DEC
               check_svd(ArmaDec,Data,prefix, suffix.c_str(),timing);
            #endif
            #ifdef _ATLAS_DEC
               check_svd(AtlasDec,Data,prefix, suffix.c_str(),timing);
            #endif
            check_svd(SapDec,Data,prefix, suffix.c_str(),timing);
      } else {
          std::string suffix=std::string("svd");
          dblarray *testRectArray =generate_random_array(5000,128,0,0,1);
          sprintf(tempname,"%s_input_SVD.fits",prefix);
          fits_write_dblarr(tempname,* testRectArray);
          std::cout << "\t Test time for SVD of matrix: [" <<
             testRectArray->nx()<<","<<testRectArray->ny()<<"]"<<std::endl;
          std::cout << "Input: "<<tempname <<std::endl;
          #ifdef _ARMA_DEC
          check_svd(ArmaDec,*testRectArray,prefix, suffix.c_str(),timing);
          #endif
          #ifdef _ATLAS_DEC
          check_svd(AtlasDec,*testRectArray,prefix, suffix.c_str(),timing);
          #endif
          check_svd(SapDec,*testRectArray,prefix, suffix.c_str(),timing);
         delete testRectArray;
      }
      //TEST SQRT MAT VIA SVD
      if(testInput){
            std::string suffix=std::string("sqrtSVD");
            #ifdef _ARMA_DEC
                check_sqrtSVD(ArmaDec,Data,prefix, suffix.c_str(),timing);
            #endif
            #ifdef _ATLAS_DEC
                check_sqrtSVD(AtlasDec,Data,prefix, suffix.c_str(),timing);
            #endif
            check_sqrtSVD(SapDec,Data,prefix, suffix.c_str(),timing);
      } else {
        dblarray* testRectArray =init_square_array_svd();
        dblarray* SqrtMat;
        #ifdef _ARMA_DEC
            std::cout << "Test Decomposition "<<ArmaDec->getSignature()<<std::endl;
            SqrtMat=ArmaDec->sqrtSVD(*testRectArray,0);
            check_sqrtSVD_simple_array(*SqrtMat);
            delete SqrtMat;
        #endif
        #ifdef _ATLAS_DEC
            std::cout << "Test Decomposition "<<AtlasDec->getSignature()<<std::endl;
            SqrtMat=AtlasDec->sqrtSVD(*testRectArray,0);
            check_sqrtSVD_simple_array(*SqrtMat);
            delete SqrtMat;
        #endif
        std::cout << "Test Decomposition "<<SapDec->getSignature()<<std::endl;
        SqrtMat=SapDec->sqrtSVD(*testRectArray,0);
        check_sqrtSVD_simple_array(*SqrtMat);
        delete SqrtMat;
      }
   }

   if(testPInvert){
      dblarray *testRectArray= init_simple_rect_array();
      std::cout << "Test Pinverse "<<SapDec->getSignature()<<std::endl;
      MatrixPInvert=SapDec->pinvert(*testRectArray,0.);
      check_pinverse_simple_rect_array(* MatrixPInvert);
      delete MatrixPInvert;
      #ifdef _ATLAS_DEC
         std::cout<< "Test Pinverse "<<AtlasDec->getSignature()<<std::endl;
         MatrixPInvert = AtlasDec->pinvert(*testRectArray);
         check_pinverse_simple_rect_array(* MatrixPInvert);
         delete MatrixPInvert;
      #endif
      #ifdef _ARMA_DEC
         std::cout << "Test Pinverse "<<ArmaDec->getSignature()<<std::endl;
         MatrixPInvert = ArmaDec->pinvert(*testRectArray);
         check_pinverse_simple_rect_array(* MatrixPInvert);
         delete MatrixPInvert;
      #endif
      delete testRectArray;
      if(testInput){
         if(Nx<Ny)
            std::cout << "Cannot compute pseudo-inverse for matrix: [" <<
                                                  Nx <<","<<Ny<<"]"<<std::endl;
         else{
            std::string suffix=std::string("pinvert");
            #ifdef _ARMA_DEC
               check_pinvert(ArmaDec,Data,prefix, suffix.c_str(),timing);
            #endif
            #ifdef _ATLAS_DEC
               check_pinvert(AtlasDec,Data,prefix, suffix.c_str(),timing);
            #endif
            check_pinvert(SapDec,Data,prefix, suffix.c_str(),timing);
         }
      } else {
          std::string suffix=std::string("pinvert");
          dblarray *testRectArray =generate_random_array(5000,128,0,0,1);
          sprintf(tempname,"%s_input_pinverse.fits",prefix);
          fits_write_dblarr(tempname,* testRectArray);
          std::cout << "\t Test time for  pseudo inverse of matrix: [" <<
             testRectArray->nx()<<","<<testRectArray->ny()<<"]"<<std::endl;
          std::cout << "Input: "<<tempname <<std::endl;
          #ifdef _ARMA_DEC
          check_pinvert(ArmaDec,* testRectArray,prefix, suffix.c_str(),timing);
          #endif
          #ifdef _ATLAS_DEC
          check_pinvert(AtlasDec,* testRectArray,prefix, suffix.c_str(),timing);
          #endif
          check_pinvert(SapDec,* testRectArray,prefix, suffix.c_str(),timing);
         delete testRectArray;
      }
   }

   delete SapDec;
   #ifdef _ATLAS_DEC
      delete AtlasDec;
   #endif
   #ifdef _ARMA_DEC
      delete ArmaDec;
      delete ArmaDecStd;
   #endif
   exit(0);
}

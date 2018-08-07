/*******************************************************************************
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
 **    File:  MatrixDecomposition.cpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Implementation of wrapper classes for several (matrix)
 **         decompositions based on different libraries.
 **    -----------
 **  - 06/2016: Sapdev, Atlas and Armadillo decompositions implemented (svd,
         inversion,pseudo inverse, cholesky except for sapdev)
 **  - 07/2016: Helper routines implemented in MatrixDecomposition
 ******************************************************************************/


#include "MatrixDecomposition.h"

#ifdef ENABLE_VERBOSE_INFO_IN_MatrixDecomposition
#define FINFO_ DEBUG_(std::cout << "> " << __PRETTY_FUNCTION__ << std::endl;)
#define FOFNI_ DEBUG_(std::cout << "< END " << __PRETTY_FUNCTION__ << std::endl;)
#else
#undef FINFO_
#undef FOFNI_
#define FINFO_
#define FOFNI_
#endif

#define DEBUGME_(x) x
//#define DEBUGME_(x)

/*******************************************************************************
 ** General MatrixDecomposition routines
*******************************************************************************/

//Start with basic operations on Columns
dblarray* MatrixDecomposition::getMeanCol2DData(dblarray &inputData){
    FINFO_;
    //InputData is row major order, we want the mean of each column
    //among the inputData.nx() columns
    size_t Ncols= inputData.nx(), Nrows=inputData.ny();
    dblarray *MeanArray= new dblarray(Ncols);
    double *meanptr=MeanArray->buffer();
    double *bufInput= inputData.buffer();
    MeanArray->init(0.);
    for(size_t krow=0; krow <Nrows;++ krow){
       meanptr=MeanArray->buffer();
       for(size_t kcol=0; kcol <Ncols;++kcol,++bufInput,++meanptr){
          *meanptr +=*bufInput;
       }
    }
    meanptr= MeanArray->buffer();
    for(size_t kcol=0; kcol <Ncols;++kcol,++meanptr){
       *meanptr /= Nrows;
    }
    FOFNI_;
    return MeanArray;
}



dblarray* MatrixDecomposition::getL2NormCol2DData(dblarray &inputData){
    FINFO_;
    //InputData is row major order, we want the l2 norm of each column
    //among the inputData.nx() columns
    size_t Ncols= inputData.nx(), Nrows=inputData.ny();
    dblarray *L2Array= new dblarray(Ncols);
    double *L2ptr= L2Array->buffer();
    double *bufInput= inputData.buffer();
    L2Array->init(0.);
    for(size_t krow=0; krow <Nrows;++ krow){
        L2ptr = L2Array->buffer();
        for(size_t kcol=0; kcol <Ncols;++kcol,++bufInput,++L2ptr){
            *L2ptr += (*bufInput)*(*bufInput);
        }
    }
    L2ptr = L2Array->buffer();
    for(size_t kcol=0; kcol <Ncols;++kcol,++L2ptr) *L2ptr = sqrt(*L2ptr);
    FOFNI_;
    return L2Array;
}

void MatrixDecomposition::addCol2DData(dblarray &inputData,dblarray &addData){
    FINFO_;
    //InputData/subData are row major order, we want to add addData from
    //each column among the inputData.ny() columns
    size_t Ncols=inputData.nx(),Nrows=inputData.ny();
    if((size_t) addData.nx() !=Ncols) {
        std::cout << "MatrixDecomposition::addCol2DData incompatible size:"<<
          " Ncol1="<< Ncols << " Ncol2="<<addData.nx ()<< std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    double *addptr= addData.buffer();
    double *bufInput= inputData.buffer();
    //Dictionary is row major order, with rows= pix, columns=Atoms
    for(size_t krow=0;krow<Nrows;++krow){
        addptr= addData.buffer();
        for(size_t kcol=0;kcol<Ncols;++kcol,++bufInput,++addptr)
                                                 *bufInput+=*addptr;
    }
    FOFNI_;
}

void MatrixDecomposition::subCol2DData(dblarray &inputData,dblarray &subData){
    FINFO_;
    //InputData/subData are row major order, we want to remove subData from
    //each column among the inputData.nx() columns
    size_t Ncols=inputData.nx(),Nrows=inputData.ny();
    if((size_t) subData.nx() !=Ncols) {
        std::cout << "MatrixDecomposition::subCol2DData incompatible size:"<<
           " Ncol1="<< Ncols << " Ncol2="<<subData.nx ()<< std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    double *subptr= subData.buffer();
    double *bufInput= inputData.buffer();
    //Dictionary is row major order, with rows= pix, columns=Atoms
    for(size_t krow=0;krow<Nrows;++krow){
        subptr= subData.buffer();
        for(size_t kcol=0;kcol<Ncols;++kcol,++bufInput,++subptr)
                                                *bufInput-=*subptr;
    }
    FOFNI_;
}

void MatrixDecomposition::multCol2DData(dblarray &inputData,dblarray &multData){
    FINFO_;
    //InputData/subData are row major order, we want to multiply by multData
    //each column among the inputData.nx() columns
    size_t Ncols=inputData.nx(),Nrows=inputData.ny();
    if((size_t) multData.nx() !=Ncols) {
        std::cout << "MatrixDecomposition::multCol2DData incompatible size:"<<
          " Ncol1="<< Ncols << " Ncol2="<<multData.nx ()<< std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    double *multptr= multData.buffer();
    double *bufInput= inputData.buffer();
    //Dictionary is row major order, with rows= pix, columns=Atoms
    for(size_t krow=0;krow<Nrows;++krow){
        multptr = multData.buffer();
        for(size_t kcol=0;kcol<Ncols;++kcol,++bufInput,++multptr)
                                                *bufInput*=*multptr;
    }
    FOFNI_;
}

void MatrixDecomposition::divCol2DData(dblarray &inputData,dblarray &divData){
    FINFO_;
    //InputData/subData are row major order, we want to divide by divData
    //each column among the inputData.nx() columns
    size_t Ncols=inputData.nx(),Nrows=inputData.ny();
    if((size_t) divData.nx() !=Ncols) {
        std::cout << "MatrixDecomposition::divCol2DData incompatible size:"<<
           " Ncol1="<< Ncols << " Ncol2="<<divData.nx ()<< std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    double *divptr= divData.buffer();
    double *bufInput= inputData.buffer();
    //Dictionary is row major order, with rows= pix, columns=Atoms
    for(size_t krow=0;krow<Nrows;++krow){
        divptr = divData.buffer();
        for(size_t kcol=0;kcol<Ncols;++kcol,++bufInput,++divptr)
                                                 *bufInput/=* divptr;
    }
    FOFNI_;
}


void MatrixDecomposition::centerCol2DData(dblarray &inputData,char *MeanName){
    FINFO_;
    //InputData is row major order, we want to center each column
    //among the inputData.nx() columns
    //Example: Dictionary in DL: where rows=one pixel, column=one atom, and
    //         we center the atoms.
    //Example: Training set in DL: where rows=one patch, column=one pixel, and
    //         we remove the mean patch
    dblarray* MeanArray=getMeanCol2DData(inputData);
    if(MeanName!=NULL) fits_write_dblarr(MeanName, *MeanArray);
    subCol2DData(inputData,*MeanArray);
    delete MeanArray;
    FOFNI_;
}

void MatrixDecomposition::l2ProjCol2DData(dblarray &inputData){
    FINFO_;
    //InputData is row major order, we want to project onto l2 ball each column
    //among the inputData.nx() columns
    //Example: Dictionary in DL: where rows=one pixel, column=one atom, and
    //         we center the atoms.
    size_t Ncols= inputData.nx(), Nrows=inputData.ny();
    dblarray *L2Array=getL2NormCol2DData(inputData);
    double *bufInput= inputData.buffer(),*L2ptr= L2Array->buffer(),*startptr;
    bufInput =inputData.buffer();
    for(size_t krow=0;krow<Nrows;++krow){
        startptr= L2ptr;
        for(size_t kcol=0;kcol<Ncols;++kcol,++bufInput,++startptr){
            if(*startptr >1) *bufInput/=(*startptr);
        }
    }
    delete L2Array;
    FOFNI_;
}

//! @remarks
//!     AW : This name is misleading, as the vectors are not normalized (l2 norm = 1), but centered and projected on the l2 unity ball
//!     Maybe rename it as l2ProjCenterCol2DData or something similar without 'norm' in the name
void MatrixDecomposition::l2NormalizeCol2DData(dblarray &inputData){
    FINFO_;
    //InputData is row major order, we want to normalize each column
    //among the inputData.nx() columns - i.e. center it and project it on l2ball
    //Example: Dictionary in DL: where rows=one pixel, column=one atom, and
    //         we center the atoms.
    centerCol2DData(inputData,NULL);
    l2ProjCol2DData(inputData);
    FOFNI_;
}

//Same Operations but with respect to rows
dblarray* MatrixDecomposition::getMeanRow2DData(dblarray &inputData){
    FINFO_;
    //InputData is row major order, we want the mean of each row
    //among the inputData.ny() rows
    size_t Ncols= inputData.nx(), Nrows=inputData.ny();
    dblarray *MeanArray= new dblarray(Nrows);
    double *meanptr=MeanArray->buffer();
    double *bufInput= inputData.buffer();
    for(size_t krow=0;krow<Nrows;++krow,++meanptr){
        *meanptr =0.;
        for(size_t kcol=0; kcol<Ncols;++kcol,++bufInput) *meanptr +=*bufInput;
        *meanptr/=Ncols;
    }
    FOFNI_;
    return MeanArray;
}

dblarray* MatrixDecomposition::getL2NormRow2DData(dblarray &inputData){
    FINFO_;
    //InputData is row major order, we want the L2 norm of each row
    //among the inputData.ny() rows
    size_t Ncols= inputData.nx(), Nrows=inputData.ny();
    dblarray *L2Array= new dblarray(Nrows);
    double *L2ptr= L2Array->buffer();
    double *bufInput= inputData.buffer();
    for(size_t krow=0;krow<Nrows;++krow,++L2ptr){
        *L2ptr=0;
        for(size_t kcol=0; kcol <Ncols;++kcol,++bufInput){
            *L2ptr += (*bufInput)*(*bufInput);
        }
        *L2ptr=sqrt(*L2ptr);
    }
    FOFNI_;
    return L2Array;
}

void MatrixDecomposition::addRow2DData(dblarray &inputData,dblarray &addData){
    FINFO_;
    //InputData is row major order, we want to add for each row addData[row]
    size_t Ncols=inputData.nx(),Nrows=inputData.ny();
    if((size_t) addData.nx() !=Nrows) {
        std::cout << "MatrixDecomposition::addRow2DData incompatible size:"<<
         " Ncol1="<< Ncols << " Ncol2="<<addData.nx()<< std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    double *addptr= addData.buffer();
    double *bufInput= inputData.buffer();
    //Dictionary is row major order, with rows= pix, columns=Atoms
    for(size_t krow=0;krow<Nrows;++krow,++ addptr){
        for(size_t kcol=0;kcol<Ncols;++kcol,++bufInput) *bufInput+=* addptr;
    }
    FOFNI_;
}

void MatrixDecomposition::subRow2DData(dblarray &inputData,dblarray &subData){
    //FINFO_;
    //InputData is row major order, we want to subtract for each row subData[row]
    size_t Ncols=inputData.nx(),Nrows=inputData.ny();
    if((size_t) subData.nx() !=Nrows) {
        std::cout << "MatrixDecomposition::subRow2DData incompatible size:"<<
         " Ncol1="<< Ncols << " Ncol2="<<subData.nx ()<< std::endl;
        MyExit::exit(EXIT_FAILURE);   }
    double *subptr= subData.buffer();
    double *bufInput= inputData.buffer();
    //Dictionary is row major order, with rows= pix, columns=Atoms
    for(size_t krow=0;krow<Nrows;++krow,++subptr){
        for(size_t kcol=0;kcol<Ncols;++kcol,++bufInput) *bufInput-=* subptr;
    }
    //FOFNI_;
}

void MatrixDecomposition::multRow2DData(dblarray &inputData,dblarray &multData){
    FINFO_;
    //InputData is row major order, we want to multiply each row by multData[row]
    size_t Ncols=inputData.nx(),Nrows=inputData.ny();
    if((size_t) multData.nx() !=Nrows) {
        std::cout << "MatrixDecomposition:: multRow2DData incompatible size:"<<
          " Ncol1="<< Ncols << " Ncol2="<<multData.nx ()<< std::endl;
        MyExit::exit(EXIT_FAILURE);   }
    double *multptr= multData.buffer();
    double *bufInput= inputData.buffer();
    //Dictionary is row major order, with rows= pix, columns=Atoms
    for(size_t krow=0;krow<Nrows;++krow,++multptr){
        for(size_t kcol=0;kcol<Ncols;++kcol,++bufInput) *bufInput*=* multptr;
    }
    FOFNI_;
}

void MatrixDecomposition::divRow2DData(dblarray &inputData,dblarray &divData){
    FINFO_;
    //InputData is row major order, we want to divide each row by divData[row]
    size_t Ncols=inputData.nx(),Nrows=inputData.ny();
    if((size_t) divData.nx() !=Nrows) {
        std::cout << "MatrixDecomposition::divRow2DData incompatible size:"<<
         " Ncol1="<< Ncols << " Ncol2="<<divData.nx ()<< std::endl;
        MyExit::exit(EXIT_FAILURE);   }
    double *divptr= divData.buffer();
    double *bufInput= inputData.buffer();
    //Dictionary is row major order, with rows= pix, columns=Atoms
    for(size_t krow=0;krow<Nrows;++krow,++ divptr){
        for(size_t kcol=0;kcol<Ncols;++kcol,++bufInput) *bufInput/=* divptr;
    }
    FOFNI_;
}

void MatrixDecomposition::centerRow2DData(dblarray &inputData){
    FINFO_;
    //InputData is row major order, we want to center each row
    //among the inputData.ny() rows
    //Example: Training set in DL: where rows=one patch, column=one pixel, and
    //         we center the patches.
    dblarray* meanArray= getMeanRow2DData(inputData);
    subRow2DData(inputData,*meanArray);
    FOFNI_;
}

void MatrixDecomposition::centerRow2DData(dblarray &inputData,
                                                          dblarray *&meanArray){
    FINFO_;
    //InputData is row major order, we want to center each row
    //among the inputData.ny() rows
    //Example: Training set in DL: where rows=one patch, column=one pixel, and
    //         we center the patches.
    meanArray= getMeanRow2DData(inputData);
    subRow2DData(inputData,*meanArray);
    FOFNI_;
}

void MatrixDecomposition::centerRow2DData(dblarray &inputData,char *MeanName){
    FINFO_;
    //InputData is row major order, we want to center each row
    //among the inputData.ny() rows
    //Example: Training set in DL: where rows=one patch, column=one pixel, and
    //         we center the patches.
    dblarray* MeanArray=getMeanRow2DData(inputData);
    if(MeanName!=NULL) fits_write_dblarr(MeanName, *MeanArray);
    subRow2DData(inputData,*MeanArray);
    delete MeanArray;
    FOFNI_;
}


//! Gets the mean along the dim-th dimension
/*!
 * InputData is row major order : Data(x,y,z) = data[z][y][x]
 * dim=0 : along 'x' : mean of each row. Output (nz,ny)
 * dim=1 : along 'y' : mean of each column. Output (nx,nz)
 * dim=2 : along 'z' : mean along third dimension. Output (nx,ny)
 */
void MatrixDecomposition::getMeanAlongDim(const dblarray &inputData,
                                          size_t dim, dblarray *meanArray) {
    FINFO_; DEBUG_(std::cout << " along dim = " << dim << std::endl;);
    inputData.printInfo("io_data");
    size_t nx = max(inputData.nx(), 1);
    size_t ny = max(inputData.ny(), 1);
    size_t nz = max(inputData.nz(), 1);
    dblarray MeanArray;
    size_t n0, n1;
    switch(dim)
    {
        case 0: // over x
            n0 = nz;
            n1 = ny;
            // => address : + kz + ky*nz;
            break;
        case 1: // over y
            n0 = nx;
            n1 = nz;
            // => address : + kx + kz*nx;
            break;
        case 2: // over z
            n0 = nx;
            n1 = ny;
            // => address : + kx + ky*nx;
            break;
    }
    if(meanArray != NULL) {
        meanArray->alloc(n0, n1, 0);
        meanArray->printInfo("*meanArray");
        MeanArray.alloc(meanArray->buffer(), n0, n1, 0);
    } else {
        MeanArray.alloc(n0, n1, 0);
    }
    DEBUGME_(MeanArray.printInfo("MeanArray"););

    double *ptrMean;
    double *ptrData = inputData.buffer();
    switch(dim) {
        case 0: // mean along x. Reading in inputData order
            for(size_t kz=0 ; kz<nz ; ++kz) {
                for(size_t ky=0 ; ky<ny ; ++ky) {
                    ptrMean = MeanArray.buffer() + kz + ky*nz;
                    *ptrMean = 0.;
                    // Evaluate Sum
                    for(size_t kx=0; kx<nx; ++kx, ++ptrData)
                        *ptrMean += *ptrData;
                    // Sum -> Mean
                    *ptrMean /= nx;
                }
            }
            break;
        case 1: // mean along y. Reading in inputData order
            MeanArray.init(0.);
            double *ptrDataKz;
            for(size_t kz=0 ; kz<nz ; ++kz) {
                // Evaluate sum
                ptrDataKz = MeanArray.buffer() + kz*nx; // + kx + kz*nx
                for(size_t ky=0 ; ky<ny ; ++ky) {
                    ptrMean = ptrDataKz;
                    for(size_t kx=0; kx<nx; ++kx, ++ptrData, ++ptrMean)
                        *ptrMean += *ptrData;
                }
                // Sum -> Mean
                ptrMean = ptrDataKz;
                for(size_t kx=0; kx<nx; ++kx, ++ptrData, ++ptrMean)
                    *ptrMean /= ny;
            }
            break;
        case 2: // mean along z. Reading in inputData order
            MeanArray.init(0.);
            // Evaluate sum
            for(size_t kz=0 ; kz<nz ; ++kz) {
                for(size_t ky=0 ; ky<ny ; ++ky) {
                    ptrMean = MeanArray.buffer() + ky*nx;// + kx;
                    for(size_t kx=0; kx<nx; ++kx, ++ptrData, ++ptrMean)
                        *ptrMean += *ptrData;
                }
            }
            // Sum -> Mean
            ptrMean = MeanArray.buffer();
            for(size_t kxy=0 ; kxy<ny*nx ; ++kxy)
                *ptrMean /= nz;
            break;
    }
    FOFNI_;
}
void MatrixDecomposition::addAlongDim(dblarray &io_data,
                                      size_t dim, dblarray MeanArray, double factor) {
    FINFO_; DEBUG_(std::cout << " along dim = " << dim << std::endl;);
    DEBUGME_(io_data.printInfo("io_data"););
    DEBUGME_(MeanArray.printInfo("MeanArray"););
    size_t nx = max(io_data.nx(), 1);
    size_t ny = max(io_data.ny(), 1);
    size_t nz = max(io_data.nz(), 1);
    double *ptrMean;
    double *ptrData = io_data.buffer();
    switch(dim) {
        case 0: // mean along x. Reading in inputData order
            double *ptrDataKzKy;
            for(size_t kz=0 ; kz<nz ; ++kz) {
                for(size_t ky=0 ; ky<ny ; ++ky) {
                    ptrMean = MeanArray.buffer() + kz + ky*nz;
                    ptrDataKzKy = ptrData;
                    // Subtract Mean
                    ptrData = ptrDataKzKy;
                    for(size_t kx=0; kx<nx; ++kx, ++ptrData)
                        *ptrData += *ptrMean*factor;
                }
            }
            break;
        case 1: // mean along y. Reading in inputData order
            MeanArray.init(0.);
            double *ptrDataKz;
            for(size_t kz=0 ; kz<nz ; ++kz) {
                ptrDataKz = MeanArray.buffer() + kz*nx; // + kx + kz*nx
                // Subtract Mean
                for(size_t ky=0 ; ky<ny ; ++ky) {
                    ptrMean = ptrDataKz;
                    for(size_t kx=0; kx<nx; ++kx, ++ptrData, ++ptrMean)
                        *ptrData += *ptrMean*factor;
                }
            }
            break;
        case 2: // mean along z. Reading in inputData order
            // Subtract Mean
            for(size_t kz=0 ; kz<nz ; ++kz) {
                for(size_t ky=0 ; ky<ny ; ++ky) {
                    ptrMean = MeanArray.buffer() + ky*nx;// + kx;
                    for(size_t kx=0; kx<nx; ++kx, ++ptrData, ++ptrMean)
                        *ptrData += *ptrMean*factor;
                }
            }
            break;
    }
    FOFNI_;
}
void MatrixDecomposition::subtractMeanAlongDim(dblarray &io_data,
                                               size_t dim, dblarray *meanArray){
    FINFO_;
    bool clearMeanArray = false;
    if(meanArray==NULL) {
        meanArray = new dblarray();
        clearMeanArray = true;
    }
    getMeanAlongDim(io_data, dim, meanArray);
    addAlongDim(io_data, dim, *meanArray, -1);
    if(clearMeanArray)
        delete meanArray;
    FOFNI_;
}

void MatrixDecomposition::l2ProjRow2DData(dblarray &inputData){
    //InputData is row major order, we want to project onto l2 ball each row
    //among the inputData.ny() rows
    //Example: Training set in DL: where rows=one patch, column=one pixel, and
    //         we center the patches.

    FINFO_;
    size_t Ncols=inputData.nx(),Nrows=inputData.ny();
    double *bufInput= inputData.buffer();
    dblarray* L2Array=getL2NormRow2DData(inputData);
    double* L2ptr= L2Array->buffer();
    //Dictionary is row major order, with rows= pix, columns=Atoms
    for(size_t krow=0; krow <Nrows;++krow,++L2ptr){
        for(size_t kcol=0; kcol<Ncols;++kcol,++bufInput) {
            if(*L2ptr>1) *bufInput/= *L2ptr;
        }
    }
    delete L2Array;
    FOFNI_;
}


void MatrixDecomposition::l2NormalizeRow2DData(dblarray &inputData){
    //InputData is row major order, we want to norm each row
    //among the inputData.ny() rows
    //Example: Training set in DL: where rows=one patch, column=one pixel, and
    //         we center the patches.
    FINFO_;
    centerRow2DData(inputData);
    l2ProjRow2DData(inputData);
    FOFNI_;
}

dblarray* MatrixDecomposition::gram(dblarray &mat){
    FINFO_;
   //Multiply tranpose(mat)##mat with only upper triangular coeff non 0
    double *matBuffer=mat.buffer();//initial position ptr for mat2
    double *itermat1,*itermat2;//initial position ptr for mat2
    double *iterGram;//initial position ptr for mat2
    dblarray *matGram;
    size_t Nx=(size_t) mat.nx(), Ny=(size_t) mat.ny();
    matGram=new dblarray(Nx,Nx,0);
    matGram->init(0.);
    for(size_t ki=0;ki<Nx;++ki){
        iterGram=matGram->buffer()+ki*Nx+ki;
        for(size_t kk=ki;kk<Nx;++kk,++iterGram) {
            itermat1= matBuffer+ki;
            itermat2= matBuffer+kk;
            for(size_t kj=0;kj<Ny;++kj,itermat1+=Nx,itermat2+=Nx)
                *iterGram +=(*itermat2)*(* itermat1);
        }
    }
    FOFNI_;
    return matGram;
}
dblarray* MatrixDecomposition::unitL2Gram(dblarray &mat){
    //Multiply tranpose(mat)##mat with only upper triangular coeff non 0
    //and mat is column l2 normalized
    double *matBuffer=mat.buffer();//initial position ptr for mat2
    double *itermat1,*itermat2;//initial position ptr for mat2
    double *iterGram;//initial position ptr for mat2
    double normcol1=0,normcol2=0;
    dblarray *matGram;
    size_t Nx=(size_t) mat.nx(), Ny=(size_t) mat.ny();
    matGram=new dblarray(Nx,Nx,0);
    matGram->init(0.);
    for(size_t ki=0;ki<Nx;++ki){
        iterGram=matGram->buffer()+ki*Nx+ki;
        for(size_t kk=ki;kk<Nx;++kk,++iterGram) {
            itermat1= matBuffer+ki;
            itermat2= matBuffer+kk;
            for(size_t kj=0;kj<Ny;++kj,itermat1+=Nx,itermat2+=Nx){
                *iterGram +=(*itermat2)*(* itermat1);
                normcol1=(*itermat1)*(*itermat1);
                normcol2=(*itermat2)*(*itermat2);
            }
            if(normcol1>0) *iterGram/=sqrt(normcol1);
            if(normcol2>0) *iterGram/=sqrt(normcol2);
        }
    }
    return matGram;
}



void  MatrixDecomposition::matMult(dblarray &Mout, dblarray &mat1, dblarray &mat2){
   //Multiply mat1##mat2, both row major
   //matProd[kj+ki*Nx2]=sum_kk mat1[kk+ki*Nx1]*mat2[kj+kk*Nx2]
    FINFO_;
   double *mat2Buffer=mat2.buffer();//initial position ptr for mat2
   double *itermat1=mat1.buffer();//ptr to iterate for mat1
   double *itermat2= mat2Buffer;//ptr to iterate for mat2
   size_t Nx1=(size_t) mat1.nx(), Ny1=(size_t) mat1.ny();
   size_t Nx2=(size_t) mat2.nx(), Ny2=(size_t) mat2.ny();
   if(Nx1!=Ny2){
      std::cout << "MatrixDecomposition:: matMult incompatible size:"<<
         " Ncol1="<< Nx1<< " Nrows2="<<Ny2<< std::endl;
      MyExit::exit(EXIT_FAILURE);
   }
   // if matProd is defined, uses the same object. otherwise object is created.
   if(((size_t)Mout.nx() != Nx2) || ((size_t)Mout.ny() != Ny1))
                                                    Mout.reform(Nx2, Ny1, 0);
   Mout.init(0.);
   double *matProdBuffer = Mout.buffer();
   double *iterprod = matProdBuffer;
   for(size_t ki=0;ki<Ny1;++ki){
      itermat2=mat2Buffer;
      for(size_t kk=0;kk<Nx1;++kk,++itermat1) {
         iterprod=matProdBuffer+ki*Nx2;
         for(size_t kj=0;kj<Nx2;++kj,++itermat2,++iterprod)
            *iterprod +=(*itermat2)*(* itermat1);
      }
   }
   FOFNI_;
}

void  MatrixDecomposition::matMult(dblarray* &Mout, dblarray &mat1, dblarray &mat2){
    if(Mout==nullptr) Mout = new dblarray();
    matMult(*Mout, mat1, mat2);
}


dblarray* MatrixDecomposition::matMult(dblarray &mat1, dblarray &mat2){
   dblarray *matProd = nullptr;
   matMult(matProd, mat1, mat2);
   return matProd;
}


void MatrixDecomposition::matMultTransp(dblarray &Mout, dblarray &mat1, dblarray &mat2) {
    //Multiply mat1##transp(mat2), both row major
    //matProd[kj+ki*Ny2]=sum_kk mat1[kk+ki*Nx1]*mat2[kk+kj*Nx2]
    double *M2Buffer = mat2.buffer(); //initial position ptr for M2
    double *itermat1 = mat1.buffer(); //ptr to iterate for M1
    double *itermat2 = M2Buffer;    //ptr to iterate for M2
    size_t Nx1=(size_t) mat1.nx(), Ny1=(size_t) mat1.ny();
    size_t Nx2=(size_t) mat2.nx(), Ny2=(size_t) mat2.ny();
    if(Nx1!=Nx2){
       std::cout << "MatrixDecomposition:: matMultTransp incompatible size:"<<
          " Ncol1="<< Nx1<< " Nrows2="<<Nx2<< std::endl;
       MyExit::exit(EXIT_FAILURE);
    }
    // if out is defined, uses the same object. otherwise  object is created.
    if(((size_t)Mout.nx() != Ny2) || ((size_t)Mout.ny() != Ny1))
                                                    Mout.reform(Ny2, Ny1, 0);
    double *matProdBuffer = Mout.buffer();
    double *iterprod = matProdBuffer;
    double *startmat1 = itermat1, *startmat2 = itermat2;
    for(size_t ki=0;ki<Ny1;++ki){
       startmat1=itermat1;
       itermat2=startmat2;
       for(size_t kj=0;kj<Ny2;++kj,++iterprod){
          *iterprod=0;
          itermat1=startmat1;
          for(size_t kk=0;kk<Nx1;++kk,++itermat1,++itermat2) {
             *iterprod +=(*itermat2)*(* itermat1);
          }
       }
    }
    FOFNI_;
}

void MatrixDecomposition::matMultTransp(dblarray* &Mout, dblarray &mat1, dblarray &mat2) {
    if(Mout==nullptr) Mout = new dblarray();
    matMultTransp(*Mout, mat1, mat2);
}

dblarray* MatrixDecomposition::matMultTransp(dblarray &mat1, dblarray &mat2){
    FINFO_;
    dblarray *matProd = nullptr;
    matMultTransp(matProd, mat1, mat2);
    FOFNI_;
    return matProd;

}

dblarray* MatrixDecomposition::transpMatMult(dblarray &mat1, dblarray &mat2){
    //Multiply mat1##mat2, both row major, and then transpose
    //matProd[ki+kj*Ny1]=sum_kk mat1[kk+ki*Nx1]*mat2[kj+kk*Nx2]
    FINFO_;
    double *mat2Buffer=mat2.buffer();//initial position ptr for mat2
    double *itermat1=mat1.buffer();//ptr to iterate for mat1
    double *itermat2= mat2Buffer;//ptr to iterate for mat2
    dblarray *matProd;
    size_t Nx1=(size_t) mat1.nx(), Ny1=(size_t) mat1.ny();
    size_t Nx2=(size_t) mat2.nx(), Ny2=(size_t) mat2.ny();
    if(Nx1!=Ny2){
        std::cout << "MatrixDecomposition:: matMultTransp incompatible size:"<<
            " Ncol1="<< Nx1<< " Nrows2="<<Ny2<< std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    matProd=new dblarray(Ny1, Nx2,0);
    matProd->init(0.);
    double *matProdBuffer=matProd->buffer();
    double *iterprod= matProdBuffer;
    for(size_t ki=0;ki<Ny1;++ki){
        itermat2=mat2Buffer;
        for(size_t kk=0;kk<Nx1;++kk,++itermat1) {
            iterprod=matProdBuffer+ki;
            for(size_t kj=0;kj<Nx2;++kj,++itermat2,iterprod+=Ny1)
            *iterprod +=(*itermat2)*(* itermat1);
        }
    }
    FOFNI_;
    return matProd;
}


dblarray* MatrixDecomposition::matRowVec(const dblarray &mat,
                                    const dblarray &vect, unsigned int index){
    //Beware, here vector is dblarray(Nelem,Nindices,0), i.e. it is composed of
    // Nindices row vector, of which the one indexed by index is used.
    //out[ki]=sum_kk mat[kk+ki*Nx1]*vect[kk+index*Nx2]
    FINFO_;
    size_t Nx1=(size_t) mat.nx(), Ny1=(size_t) mat.ny();
    size_t Nx2=(size_t) vect.nx(), Ny2=(size_t) vect.ny();

    dblarray * outRowVect;
    if(index>Ny2){
        std::cout << "MatrixDecomposition:: matRowVec not enough row vectors:"<<
                    index<<" index vs Nrows="<< Ny2<<std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    if(Nx1!=Nx2){
        std::cout << "MatrixDecomposition:: matRowVec incompatible size:"<<
            " Ncol1="<< Nx1<< " Nrows2="<<Nx2<< std::endl;
        MyExit::exit(EXIT_FAILURE);
    }
    double *itermat=mat.buffer();
    double *vectBuffer=vect.buffer()+index*Nx2;
    double *itervect= vectBuffer;

    outRowVect=new dblarray(Ny1,0,0);
    double *iterout = outRowVect->buffer();
    for(size_t ki=0;ki<Ny1;++ki,++iterout){
        *iterout=0.;
        itervect = vectBuffer;
        for(size_t kk=0;kk<Nx1;++kk,++itermat,++itervect)
            *iterout +=(*itermat)*(*itervect);
    }
    //outVect[ki]=sum_kk mat1[kk+ki*Nx1]*mat2[kk]
    FOFNI_;
    return outRowVect;
}


double MatrixDecomposition::dotProduct(const double* vector1,
            const double* vector2, size_t stride1,size_t stride2,size_t Nelems,
             const double* metric,bool useMetric,bool diagMetric){
    FINFO_;
    double dotProduct=0.;
    if(!useMetric){
        for(size_t kel1=0,kel2=0;kel1<Nelems*stride1;kel1+=stride1,kel2+=stride2)
            dotProduct+=vector1[kel1]*vector2[kel2];
    } else if (diagMetric){
        const double *metricBuffer=metric;
        for(size_t kel1=0,kel2=0;kel2<Nelems*stride2;kel1+=stride1,kel2+=stride2,
        ++metricBuffer)
                    dotProduct+= vector1[kel1] * (*metricBuffer) * vector2[kel2];
    } else {
        double transvec=0.;
        const double *metricBuffer=metric;
        const double *buff1=vector1,*buff2=vector2;
        for(size_t kvec=0;kvec<Nelems;++kvec,buff2+=stride2){
            transvec=0.;
            buff1=vector1;
            for(size_t kel=0;kel<Nelems;kel++,buff1+=stride1,++metricBuffer)
            transvec+= (*buff1) * (*metricBuffer);
            dotProduct+=transvec* (*buff2);
        }
    }
    FOFNI_;
    return dotProduct;
}


double MatrixDecomposition::dotProduct(const dblarray &vector1,
                              const dblarray &vector2, const dblarray &metric,
                                                bool useMetric,bool diagMetric){
    FINFO_;
    if(vector1.nx() == vector2.nx()){
        FOFNI_;
        return MatrixDecomposition::dotProduct(vector1.buffer(),vector2.buffer(),
            1,1,vector1.nx(),metric.buffer(),useMetric,diagMetric);
    } else {
         std::cout << "MatrixDecomposition::Input dimensions do not match: " <<
             vector1.nx() << " vs "<< vector2.nx()<< std::endl;
         MyExit::exit(EXIT_FAILURE);
         return 0;
    }
    FOFNI_;
}

void MatrixDecomposition::solve_UT(dblarray& A, dblarray& b, dblarray& x) {
    std::cout << "MatrixDecomposition::solve_UT not implemented yet. Use Armadillo."<< std::endl;
    MyExit::exit(EXIT_FAILURE);
}
void MatrixDecomposition::solve_LT(dblarray& A, dblarray& b, dblarray& x) {
    std::cout << "MatrixDecomposition::solve_LT not implemented yet. Use Armadillo."<< std::endl;
    MyExit::exit(EXIT_FAILURE);
}

dblarray* MatrixDecomposition::sqrtSVD(dblarray& Mat, double relthr) {
    FINFO_;
    dblarray U,S,Vt, *sqrtMat;
    //Mat[m_Nx,m_Nx], sqrtMat[m_N,m_Nx]
    //U[nSing,m_Ny], Vt[m_Nx,nSing], S[nSing]
    //Ut[m_Ny,nSing], V[nSing,m_Nx]
    size_t m_Nx=Mat.nx();
    size_t m_Ny=Mat.ny();
    if(m_Ny!=m_Nx){
        std::cout << "MatrixDecomposition::sqrtSVD dimensions do not match: " <<
            m_Nx << " != "<< m_Ny<< std::endl;
        MyExit::exit(EXIT_FAILURE);
    }

    svd(Mat,U,S,Vt);
    sqrtMat=new dblarray(m_Nx, m_Nx,0);
    sqrtMat->init(0.);
    double thr= relthr*fabs(S.maxfabs());
    if(thr<0) {
        printf("MatrixDecomposition::Wrong value for thr %lf, should be >0\n",thr);
        MyExit::exit(EXIT_FAILURE);
    }
    size_t NxU=U.nx(),nSing=S.nx();
    double *bufS= S.buffer();

    //sqrt singular vectors
    bufS= S.buffer();
    for(size_t ksing=0;ksing<nSing;++ksing,++bufS) {
        //std::cout<<"SINGVAL["<<ksing<<"]="<<*bufS<<std::endl;
        if(fabs(*bufS)>thr) *bufS=sqrt(fabs(*bufS));
        else *bufS=0.;
    }
    double *bufU;
    //First U<-U##sqrt(S) ; U[j+i*nSing]=s[j] U[j+i*m_Nx]
    for(size_t ki=0;ki<m_Ny;++ki){
        bufU=U.buffer()+ki* NxU;//Line Ui
        bufS= S.buffer();
        for(size_t kj=0;kj<nSing;++kj,++bufU,++bufS) *bufU*=(*bufS);
    }
    double *bufsqrt, *bufV, *startVt=Vt.buffer(),*startU=U.buffer();
    //Now multiply U##Vt, sqrtMat[kj+ki*m_Ny]=U[kk+ki*m_Ny]##Vt[kj+kk*nSing]
    for(size_t ki=0;ki<m_Nx;ki++) {
        bufV=startVt;
        for(size_t kk=0;kk<nSing;++kk) {
            bufU=startU+kk+ki*m_Ny;
            bufsqrt=sqrtMat->buffer()+ki*m_Ny;//start of line ki for pinv
            for(size_t kj=0;kj<m_Ny;++kj,++bufV,++bufsqrt)
                                                    (*bufsqrt)+=*(bufV)*(*bufU);
        }
    }
    U.free();
    Vt.free();
    S.free();
    FOFNI_;
    return sqrtMat;
}


/*******************************************************************************
 ** SapDecomposition routines
*******************************************************************************/

void SapDecomposition::svd(dblarray& Mat,dblarray& U,dblarray& S,dblarray& Vt){
    FINFO_;
    double* m_matbuffer=Mat.buffer();
    int m_naxis=Mat.naxis();
    if(m_matbuffer !=NULL) {
        if(m_naxis != 2) {
            std::cout << "SapDecomposition:: SVD only for 2dim matrices ("<<
                        Mat.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        m_matoper.svd(Mat,U,S,Vt);
    }
    FOFNI_;
}

dblarray* SapDecomposition::invert(dblarray& squareMat) {
    FINFO_;
    double* m_matbuffer=squareMat.buffer();
    size_t m_Nx=squareMat.nx();
    size_t m_Ny=squareMat.ny();
    int m_naxis=squareMat.naxis();
    dblarray* invSquareMat =new dblarray(m_Ny, m_Nx,0);
    if(m_matbuffer !=NULL) {
        if(m_naxis != 2) {
            std::cout << "SapDecomposition:: Invert only for 2dim matrices ("<<
                        squareMat.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        if(m_Nx!=m_Ny){
            std::cout<< "SapDecomposition:: Non invertible rectangular matrix ["<<
                    m_Nx <<","<<m_Ny << "]"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        m_matoper.inv_mat_svd(squareMat, *invSquareMat);
    }
    FOFNI_;
    return invSquareMat;
}

dblarray* SapDecomposition::pinvert(dblarray&  Mat, double relthr) {
    FINFO_;
    double* m_matbuffer=Mat.buffer();
    size_t m_Nx=Mat.nx();
    size_t m_Ny=Mat.ny();
    int m_naxis=Mat.naxis();

    dblarray* pinvMat =new dblarray(m_Ny, m_Nx,0);
    if(m_matbuffer !=NULL) {
        if(m_naxis != 2) {
            std::cout << "SapDecomposition::pinvert only for 2dim matrices ("<<
                        Mat.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        m_matoper.inv_mat_svd(Mat, *pinvMat, relthr);
    }
    FOFNI_;
    return pinvMat;
}

dblarray* SapDecomposition::cholesky(dblarray& PosDefMat) {
    FINFO_;
    std::cout << "Not implemented yet" <<std::endl;
    FOFNI_;
    return &(PosDefMat);
}

/*******************************************************************************
 ** AtlasDecomposition routines
*******************************************************************************/
#ifdef _ATLAS_DEC
void AtlasDecomposition::svd(dblarray& Mat, dblarray& U, dblarray& S,
                                                               dblarray& Vt) {
    FINFO_;
    double* m_matbuffer=Mat.buffer();
    size_t m_Nx=Mat.nx();
    size_t m_Ny=Mat.ny();
    int m_naxis=Mat.naxis();

    if(m_matbuffer !=NULL) {
        if(m_naxis != 2) {
            std::cout << "AtlasDecomposition:: Invert only for 2dim matrices ("<<
                        Mat.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }

        char jobz='s';
        double *tempbuffer=NULL;
        tempbuffer=(double*) malloc(m_Nx* m_Ny*sizeof(double));
        if(tempbuffer==NULL){
            std::cout<<"AtlasDecomposition::svd can't allocate memory"<<std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        size_t NSingVal= (m_Nx < m_Ny)? m_Nx : m_Ny;
        U.reform(NSingVal, m_Ny,0);
        Vt.reform(m_Nx, NSingVal,0);
        S.reform(NSingVal);
//      for(size_t kx=0,kpix=0;kx<m_Nx;++kx) {//columns
//         for(size_t offset=kx*m_Ny;offset<(kx+1)*m_Ny;++offset,++kpix) //rows
        for(size_t kpix=0;kpix<m_Nx*m_Ny;++kpix) {//columns
            tempbuffer[kpix]=m_matbuffer[kpix];
        }
        int info = 0;

        info=LAPACKE_dgesdd(LAPACK_ROW_MAJOR,jobz,m_Ny,m_Nx,tempbuffer,m_Nx,
                                S.buffer(),U.buffer(),NSingVal,Vt.buffer(),m_Nx);
        if(info) {
            std::cout<<"AtlasDecomposition::svd can't compute SVD"<<std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        free(tempbuffer);
    }
    FOFNI_;
}

dblarray* AtlasDecomposition::invert(dblarray& squareMat) {
    FINFO_;
    double* m_matbuffer=squareMat.buffer();
    size_t m_Nx=squareMat.nx();
    size_t m_Ny=squareMat.ny();
    int m_naxis=squareMat.naxis();

    if(m_matbuffer !=NULL) {
        if(m_naxis != 2) {
            std::cout << "AtlasDecomposition:: Invert only for 2dim matrices ("<<
                        squareMat.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        if(m_Nx!=m_Ny){
            std::cout<< "AtlasDecomposition:: Invert only for square matrices ["<<
                    m_Nx <<","<<m_Ny << "]"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        dblarray* invSquareMat =new dblarray(m_Nx, m_Nx,0);
        double *invA= invSquareMat->buffer();
        size_t Nx2=(size_t) m_Nx * m_Nx;
        int *IPIV = new int[m_Nx +1];
        //first copy invA
        for(size_t kpix=0;kpix<Nx2;kpix++) invA[kpix]= m_matbuffer[kpix];
        if(LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m_Nx, m_Nx,invA, m_Nx,IPIV)) {
            printf("Cannot perform LU factorization\n");
            MyExit::exit(EXIT_FAILURE);
        }
        if(LAPACKE_dgetri(LAPACK_ROW_MAJOR, m_Nx,invA, m_Nx,IPIV)){
            printf("Cannot invert matrix\n");
            MyExit::exit(EXIT_FAILURE);
        }
        delete IPIV;
        FOFNI_;
        return invSquareMat;
    }
    FOFNI_;
    return NULL;
}

dblarray* AtlasDecomposition::pinvert(dblarray& Mat, double relthr) {
    FINFO_;
    dblarray U,S,Vt, *pinv;
    //Mat[m_Nx,m_Ny], pinv[m_Ny,m_Nx]
    //U[nSing,m_Ny], Vt[m_Nx,nSing], S[nSing]
    //Ut[m_Ny,nSing], V[nSing,m_Nx]
    size_t m_Nx=Mat.nx();
    size_t m_Ny=Mat.ny();

    AtlasDecomposition::svd(Mat,U,S,Vt);
    pinv=new dblarray(m_Ny, m_Nx,0);
    pinv->init(0.);
    double thr= relthr*fabs(S.maxfabs());
    if(thr<0) {
        printf("AtlasDecomposition::Wrong value for thr %lf, should be >0\n",thr);
        MyExit::exit(EXIT_FAILURE);
    }
    size_t NxU=U.nx(),nSing=S.nx();

    //invert singular vectors
    double *bufS= S.buffer();
    for(size_t ksing=0;ksing<nSing;++ksing,++bufS) {
        if(fabs(*bufS)>thr) *bufS=1./(*bufS);
        else *bufS=0.;
    }
    double *bufU;
    //First U<-U##1/S ; U[j+i*nSing]=s[j] U[j+i*m_Nx]
    for(size_t ki=0;ki<m_Ny;++ki){
        bufU=U.buffer()+ki* NxU;//Line Uj
        bufS= S.buffer();
        for(size_t kj=0;kj<nSing;++kj,++bufU,++bufS) *bufU*=(*bufS);
    }
    double *bufpinv, *bufV, *startVt=Vt.buffer(),*startU=U.buffer();
    //Now multiply V##Ut, pinv[kj+ki*m_Ny]=V[kk+ki*nSing]*Ut[kj+kk*m_Ny]
    //                                    =Vt[ki+kk*m_Nx]*U[kk+kj*nSing];
    for(size_t ki=0;ki<m_Nx;ki++) {
        for(size_t kk=0;kk<nSing;++kk) {
            bufV=startVt+kk*m_Nx+ki;
        bufU=startU+kk;
        bufpinv=pinv->buffer()+ki*m_Ny;//start of line ki for pinv
        for(size_t kj=0;kj<m_Ny;++kj,bufU+=nSing,++bufpinv)
                                                    (*bufpinv)+=*(bufV)*(*bufU);
        }
    }
    U.free();
    Vt.free();
    S.free();
    FOFNI_;
    return pinv;
}


dblarray* AtlasDecomposition::cholesky(dblarray& PosDefMat) {
    FINFO_;
    int info=0;
    double* m_matbuffer=PosDefMat.buffer();
    size_t m_Nx=PosDefMat.nx();
    size_t m_Ny=PosDefMat.ny();
    int m_naxis=PosDefMat.naxis();

    if(m_matbuffer !=NULL) {
        if(m_naxis != 2) {
            std::cout << "AtlasDecomposition::Cholesky only for 2dim matrices ("<<
                        PosDefMat.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        if(m_Nx!=m_Ny){
            std::cout<<"AtlasDecomposition:: Cholesky only for square matrices ["<<
                    m_Nx <<","<<m_Ny << "]"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        dblarray* choleskyMat =new dblarray(m_Nx, m_Nx,0);
        double *cholA= choleskyMat->buffer();
        size_t Nx2=(size_t) m_Nx * m_Nx;
        //first copy invA
        for(size_t kpix=0;kpix<Nx2;kpix++) cholA[kpix]= m_matbuffer[kpix];
        info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR,CblasUpper, m_Nx, cholA, m_Nx);
        if(info >0) {
            printf("Leading minor %d is not positive definite\n",info);
            MyExit::exit(EXIT_FAILURE);
        }
        if(info <0) {
            printf("Argument %d has illegal value\n",info);
            MyExit::exit(EXIT_FAILURE);
        }
        FOFNI_;
        return choleskyMat;
    }
    FOFNI_;
    return NULL;
}
#endif // _ATLAS_DEC

/*******************************************************************************
 ** ArmadilloDecomposition routines
*******************************************************************************/
#ifdef _ARMA_DEC
void ArmadilloDecomposition::svd(dblarray& Mat, dblarray& U, dblarray& S,
                                                                dblarray& Vt) {
    FINFO_;
    double* m_matbuffer=Mat.buffer();
    size_t m_Nx=Mat.nx();
    size_t m_Ny=Mat.ny();
    int m_naxis=Mat.naxis();
    if(m_matbuffer !=NULL) {
        if(m_naxis != 2) {
            std::cout << "ArmaDecomposition:: Invert only for 2dim matrices ("<<
                        Mat.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        bool status =true;
        arma::mat ArmaMat=arma::mat(m_Ny, m_Nx);
        //Beware :dblarr buffer is row major, but mat is column-major
        for(size_t krow=0, kpix=0;krow<m_Ny;++krow)
            for(size_t kcol=0;kcol<m_Nx;++kcol,++kpix)
                                    ArmaMat.at(krow,kcol)= m_matbuffer[kpix];
        size_t NSingVal= (m_Nx < m_Ny)? m_Nx : m_Ny;
        S.reform(NSingVal);
        U.reform(NSingVal, m_Ny,0);
        Vt.reform(m_Nx, NSingVal,0);
        dblarray Ut=dblarray(m_Ny, NSingVal,0);

        arma::mat ArmaU=arma::mat(Ut.buffer(),m_Ny, NSingVal, false,true);
        arma::mat ArmaV=arma::mat(Vt.buffer(), NSingVal, m_Nx, false,true);
        arma::vec ArmaS=arma::vec(S.buffer(), NSingVal,false,true);
        if(m_svdMethod== STD) status=arma::svd_econ(ArmaU,ArmaS,ArmaV,ArmaMat,
                                                                    "both","std");
        else status=arma::svd_econ( ArmaU, ArmaS, ArmaV, ArmaMat,"both");

        if (!status){
            std::cout << "ArmaDecomposition:: SVD failed"<<std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        //S is already allright, V is Vt, needs to transpose U)
        m_matoper.transpose(Ut,U);

    }
    FOFNI_;
}

dblarray* ArmadilloDecomposition::cholesky(dblarray& PosDefMat) {
    FINFO_;
    double* m_matbuffer=PosDefMat.buffer();
    size_t m_Nx=PosDefMat.nx();
    size_t m_Ny=PosDefMat.ny();
    int m_naxis=PosDefMat.naxis();

    if(m_matbuffer !=NULL) {
        if(m_naxis != 2) {
            std::cout << "ArmaDecomposition::Cholesky only for 2dim matrices ("<<
                        PosDefMat.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        if(m_Nx!=m_Ny){
            std::cout<<"ArmaDecomposition:: Cholesky only for square matrices ["<<
                    m_Nx <<","<<m_Ny << "]"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        arma::mat *ArmaMat= new arma::mat(m_Nx, m_Nx);
        //Beware :dblarr buffer is row major, but mat is column-major
        for(size_t krow=0, kpix=0;krow<m_Nx;krow++)
            for(size_t kcol=0;kcol<m_Nx;kcol++,kpix++)
                                    ArmaMat->at(krow,kcol)= m_matbuffer[kpix];

        arma::mat UpperMat;
        bool status=arma::chol( UpperMat, *ArmaMat );
        if(!status){
            printf("Cannot perform cholesky decomposition\n");
            MyExit::exit(EXIT_FAILURE);
        }
        delete ArmaMat;
        dblarray* choleskyMat =new dblarray(m_Nx, m_Nx,0);
        double * bufferMat= choleskyMat->buffer();
        for(size_t krow=0, kpix=0;krow<m_Nx;krow++)
            for(size_t kcol=0;kcol<m_Nx;kcol++,kpix++)
                                    bufferMat[kpix]= UpperMat.at(krow,kcol);
        FOFNI_;
        return choleskyMat;
    }
    FOFNI_;
    return NULL;
}

dblarray* ArmadilloDecomposition::invert(dblarray& squareMat) {
    FINFO_;
    double* m_matbuffer=squareMat.buffer();
    size_t m_Nx=squareMat.nx();
    size_t m_Ny=squareMat.ny();
    int m_naxis=squareMat.naxis();

    if(m_matbuffer !=NULL) {
        if(m_naxis != 2) {
            std::cout << "ArmaDecomposition:: Invert only for 2dim matrices ("<<
                        squareMat.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        if(m_Nx!=m_Ny){
            std::cout << "ArmaDecomposition:: Invert only for square matrices ["<<
                    m_Nx <<","<<m_Ny << "]"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        dblarray* invSquareMat =new dblarray(m_Nx, m_Nx,0);
        double *invA= invSquareMat->buffer();

        arma::mat* ArmaMat=new arma::mat(m_Nx, m_Nx);
        //Beware :dblarr buffer is row major, but mat is column-major
        for(size_t krow=0, kpix=0;krow<m_Nx;krow++)
            for(size_t kcol=0;kcol<m_Nx;kcol++,kpix++)
                                    ArmaMat->at(krow,kcol)= m_matbuffer[kpix];

        arma::mat invMat;
        bool status=arma::inv( invMat, *ArmaMat );
        if(!status){
            printf("Cannot perform inversion\n");
            MyExit::exit(EXIT_FAILURE);
        }
        delete ArmaMat;
        for(size_t krow=0, kpix=0;krow<m_Nx;krow++)
            for(size_t kcol=0;kcol<m_Nx;kcol++,kpix++)
                                    invA[kpix]= invMat.at(krow,kcol);
        FOFNI_;
        return invSquareMat;
    }
    FOFNI_;
    return NULL;
}

dblarray* ArmadilloDecomposition::pinvert(dblarray& Mat,double thr) {
    FINFO_;
    double* m_matbuffer=Mat.buffer();
    size_t m_Nx=Mat.nx();
    size_t m_Ny=Mat.ny();
    int m_naxis=Mat.naxis();
    if(m_matbuffer !=NULL) {
        if(m_naxis != 2) {
            std::cout << "ArmaDecomposition:: Invert only for 2dim matrices ("<<
                        Mat.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
        dblarray* pinvMat =new dblarray(m_Ny, m_Nx,0);
        double *pinvA= pinvMat->buffer();

        arma::mat* ArmaMat=new arma::mat(m_Ny, m_Nx);
        //Beware :dblarr buffer is row major, but mat is column-major
        for(size_t krow=0, kpix=0;krow<m_Ny;++krow)
            for(size_t kcol=0;kcol<m_Nx;++kcol,++kpix)
                                    ArmaMat->at(krow,kcol)= m_matbuffer[kpix];

        arma::mat ArmapinvMat;
        bool status=arma::pinv( ArmapinvMat, *ArmaMat, thr);
        if(!status){
            char tempname[1024];
            printf("Cannot perform pseudo-inversion\n");
            std::cout<<"Matrix is "<<m_Nx<<"x"<<m_Ny<<std::endl;
            sprintf(tempname,"FailedMatrixpInvert.fits");
            fits_write_dblarr(tempname, Mat);
            MyExit::exit(EXIT_FAILURE);
        }
        delete ArmaMat;
        for(size_t krow=0, kpix=0;krow<m_Nx;++krow)
            for(size_t kcol=0;kcol<m_Ny;++kcol,++kpix)
                                    pinvA[kpix]= ArmapinvMat.at(krow,kcol);
        FOFNI_;
        return pinvMat;
    }
    FOFNI_;
    return NULL;
}

void ArmadilloDecomposition::solve_UT(dblarray& A, dblarray& b, dblarray& x) {
    double *A_p = A.buffer();
    size_t Nx = A.nx();
    size_t Ny =A.ny();
    size_t nAxis =A.naxis();
    if(A_p !=NULL) {
        if(nAxis != 2) {
            std::cout << "ArmaDecomposition:: solve only for 2dim matrices ("<<
                         A.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }

        bool status = true;

        //Beware: dblarr buffer is row major, but arma::mat is column-major
        arma::mat ArmaA=arma::mat(Ny, Nx);
        for(size_t krow=0, kpix=0;krow<Ny;++krow)
            for(size_t kcol=0;kcol<Nx;++kcol,++kpix)
                ArmaA.at(krow,kcol)= A_p[kpix];

        if((size_t) x.n_elem() != Nx)
            x.reform(Nx);
        arma::mat ArmaB = arma::vec(b.buffer(), Ny, false, true);
        arma::vec ArmaX = arma::vec(x.buffer(), Nx, false, true);

        status = arma::solve(ArmaX, arma::trimatu(ArmaA), ArmaB);

        if (!status){
            std::cout << "ArmaDecomposition:: SVD failed"<<std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
    }
}

void ArmadilloDecomposition::solve_LT(dblarray& A, dblarray& b, dblarray& x) {
    double *A_p = A.buffer();
    size_t Nx = A.nx();
    size_t Ny =A.ny();
    size_t nAxis =A.naxis();
    if(A_p !=NULL) {
        if(nAxis != 2) {
            std::cout << "ArmaDecomposition:: solve only for 2dim matrices ("<<
                         A.naxis() << ")"<< std::endl;
            MyExit::exit(EXIT_FAILURE);
        }

        bool status = true;

        //Beware: dblarr buffer is row major, but arma::mat is column-major
        arma::mat ArmaA=arma::mat(Ny, Nx);
        for(size_t krow=0, kpix=0;krow<Ny;++krow)
            for(size_t kcol=0;kcol<Nx;++kcol,++kpix)
                ArmaA.at(krow,kcol)= A_p[kpix];

        if((size_t) x.n_elem() != Nx)
            x.reform(Nx);
        arma::mat ArmaB = arma::vec(b.buffer(), Ny, false, true);
        arma::vec ArmaX = arma::vec(x.buffer(), Nx, false, true);

        status = arma::solve(ArmaX, arma::trimatl(ArmaA), ArmaB);

        if (!status){
            std::cout << "ArmaDecomposition:: SVD failed"<<std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
    }
}

#endif //_ARMA_DEC

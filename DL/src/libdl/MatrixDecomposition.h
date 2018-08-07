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
 **    File:  MatrixDecomposition.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces wrapper classes for several (matrix)
 **         decompositions based on different libraries.
 **    -----------
 **  - 06/2016: Sapdev, Atlas and Armadillo decompositions added
 **  - 07/2016: Helper routines added in MatrixDecomposition
 **  - 11/2016: New Helper routines added to help readability. Documentation.
 ******************************************************************************/


#ifndef MatrixDecomposition_H
#define MatrixDecomposition_H
#define EpsDoublePrec 2e-16 //2.220446049250313e-16 //based pm Matlab
#define EpsSinglePrec 2e-7 //for floating points
#include <cmath>
//Need to define complex beforehand, otherwise mixing of "complex.h"
//and <complex> which are incompatible and cause compiler errors
#include <complex>
#include <algorithm>
#include <stdexcept>

#if defined _ARMA_DEC || defined _ATLAS_DEC
    #define lapack_complex_float std::complex<float>
    #define lapack_complex_double std::complex<double>
//   #include "lapacke.h"
//   extern "C"  {
//      #include "cblas.h"
//      #include "clapack.h"
//   }
    #define ARMA_DONT_USE_CXX11
    #ifndef _ATLAS_DEC
        #define ARMA_USE_LAPACK
        //#define ARMA_USE_BLAS
    #else //defined _ATLAS_DEC
        #include "lapacke.h"
        extern "C"  {
            #include "cblas.h"
            #include "clapack.h"
        }
    #endif
#endif


#include "TempArray.h"
#include "MatrixOper.h"
#include "IM_IO.h"
#include "myexit.hpp"
#include "ISAP_verbose.hpp"

#ifdef _ARMA_DEC
    #include <armadillo>
#endif
enum ArmadilloSVD {DC ,STD};//Standard, Divide&Conquer (Faster,output different)

//Note: dblarray is row major ordering of the buffer if(x,y)=(col,row).
/**
 * @class MatrixDecomposition
 * @brief Abstract class essentially implementing the basic helper routines.
*/
class MatrixDecomposition
{
    protected:
//          NON THREAD SAFE
//        double* m_matbuffer;/**< Ptr to matrix buffer.*/
//        size_t m_Nx;/**< Number of columns.*/
//        size_t m_Ny;/**< Number of rows.*/
//        int m_naxis;/**< Number of dimensions.*/

        double m_invcondnb;/**< minimum inverse condition number (inversion).*/
    public:
        //MatrixDecomposition():m_matbuffer(NULL), m_Nx(0), m_Ny(0),m_naxis(0),
        //                                                        m_invcondnb(0.){}
        /**
        * Standart constructor.
        */
        MatrixDecomposition():m_invcondnb(0.){}
        /**
        * Standart destructor.
        * @return number of DL iterations.
        */
        virtual ~MatrixDecomposition(){}
        //members
        MatOper m_matoper;/**< Standard ISAP matrix operator.*/
        //setters
        /**
        * Set the current array we want to work on.
        * @param[in] Mat array where this class operates on.
        */
        //NON THREAD SAFE
        /*void setPtrMatrix(const dblarray& Mat) {
            m_matbuffer=Mat.buffer();
            m_Nx=Mat.nx();
            m_Ny=Mat.ny();
            m_naxis=Mat.naxis();
        }*/
        /**
        * Set the minimal inverse condition number for inversion.
        * @param[in] maxCondNumber maximal tolerated condition number.
        */
        void setmaxCondNumber(double maxCondNumber){m_invcondnb =1/maxCondNumber;}
        /**
        * Set the minimal inverse condition number for inversion.
        * @param[in] relepseigen minimal relative ratio in between eigenvalues.
        */
        void setEigenRelThresh(double relepseigen){m_invcondnb = relepseigen;}
        //getters
        /**
        * Interface to get the signature of the object.
        * @return signature of the object
        */
        virtual std::string getSignature() const =0;
        /**
        * Get the number of columns of the array considered.
        * @return number of columns.
        */
        //size_t getNCols() const{return m_Nx;}
        /**
        * Get the number of rows of the array considered.
        * @return number of rows.
        */
        //size_t getNRows() const{return m_Ny;}
        /**
        * Get the number of dimensions of the array considered.
        * @return number of dimensions.
        */
        //size_t getNAxis() const{return m_naxis;}
        //processing functions
        /**
        * Helper routine to get the mean of a column array.
        * @param[in] inputData row major data to compute the mean.
        * @return array of means on the columns [size Ncols].
        */
        static dblarray* getMeanCol2DData(dblarray &inputData);
        /**
        * Helper routine to get the L2 norm of a column array.
        * @param[in] inputData row major data to compute the L2 norm.
        * @return array of L2 norms on the columns [size Ncols].
        */
        static dblarray* getL2NormCol2DData(dblarray &inputData);
        /**
        * Helper routine to calculate the mean of data along a given dimension.
        * @param[in] inputData row major data to update.
        * @param[in] dimension where to estimate the mean (0:x,1:y,2:z)
        * @param[out] output array of means (size depending on dimensions)
        *               /!\ if passed, must be allocated by user
        */
        static void getMeanAlongDim(const dblarray &inputData, size_t dim, dblarray *meanArray);
        /**
        * Helper routine to add an array to each column of another.
        * @param[in/out] io_data row major data to update.
        * @param[in] addData row major data to add [size Ncols] to each column.
        * @param[in] input array of means (size depending on dimensions)
        * @param[in] factor : adding factor*mean
        */
        static void addAlongDim(dblarray &io_data, size_t dim, dblarray MeanArray, double factor = 1.);
        /**
        * Helper routine to subtract the mean of data along a given dimension.
        * @param[in/out] io_data row major data to update.
        * @param[in] dimension where to remove the mean (0:x,1:y,2:z)
        * @param[out] output array of means (size depending on dimensions), must be allocated by user
        */
        static void subtractMeanAlongDim(dblarray &io_data, size_t dim, dblarray *meanArray = nullptr);
        /**
        * Helper routine to add an array to each column of another.
        * @param[in] inputData row major data to update.
        * @param[in] addData row major data to add [size Ncols] to each column.
        */
        static void addCol2DData(dblarray &inputData,dblarray &addData);
        /**
        * Helper routine to sutract an array to each column of another.
        * @param[in] inputData row major data to update.
        * @param[in] subData row major data to subtract [size Ncols]
            to each column.
        */
        static void subCol2DData(dblarray &inputData,dblarray &subData);
        /**
        * Helper routine to multiply an array to each column of another.
        * @param[in] inputData row major data to update.
        * @param[in] multData row major data to multiply [size Ncols]
                to each column.
        */
        static void multCol2DData(dblarray &inputData,dblarray &multData);
        /**
        * Helper routine to divide an array to each column of another.
        * @param[in] inputData row major data to update.
        * @param[in] divData row major data to divide [size Ncols]
                to each column.
        */
        static void divCol2DData(dblarray &inputData,dblarray &divData);
        /**
        * Helper routine to center an array (mean for each column=0).
        * @param[in] inputData row major data to update.
        * @param[in] MeanName string with filename to save fits on.
        */
        static void centerCol2DData(dblarray &inputData,char *MeanName=NULL);
        /**
        * Helper routine to project an array on l2 ball (column l2norm <1).
        * @param[in] inputData row major data to update.
        */
        static void l2ProjCol2DData(dblarray &inputData);

        /**
        * Helper routine to center and  project an array on l2 ball
            (column l2norm <1).
        * @param[in] inputData row major data to update.
        */
        static void l2NormalizeCol2DData(dblarray &inputData);
        /**
        * Helper routine to get the mean of a row array.
        * @param[in] inputData row major data to compute the mean.
        * @return array of means on the rows [size Nrows].
        */
        static dblarray* getMeanRow2DData(dblarray &inputData);
        /**
        * Helper routine to get the L2 norm of a row array.
        * @param[in] inputData row major data to compute the L2 norm.
        * @return array of L2 norms on the rows [size Nrows].
        */
        static dblarray* getL2NormRow2DData(dblarray &inputData);
        /**
        * Helper routine to add an array to each row of another.
        * @param[in] inputData row major data to update.
        * @param[in] addData row major data to add [size Nrows] to each row.
        */
        static void addRow2DData(dblarray &inputData,dblarray &addData);
        /**
        * Helper routine to subtract an array to each row of another.
        * @param[in] inputData row major data to update.
        * @param[in] subData row major data to subtract [size Nrows]
            to each row.
        */
        static void subRow2DData(dblarray &inputData,dblarray &subData);
        /**
        * Helper routine to multiply an array to each row of another.
        * @param[in] inputData row major data to update.
        * @param[in] multData row major data to multiply [size Nrows]
            to each row.
        */
        static void multRow2DData(dblarray &inputData,dblarray &multData);
        /**
        * Helper routine to divide an array to each row of another.
        * @param[in] inputData row major data to update.
        * @param[in] divData row major data to divide [size Nrows]
            to each row.
        */
        static void divRow2DData(dblarray &inputData,dblarray &divData);
        /**
        * Helper routine to center an array (mean for each row=0).
        * @param[in] inputData row major data to update.
        * @param[out] meanArray array with means for each row [size Nrows].
        */
        static void centerRow2DData(dblarray &inputData,
                                                    dblarray *&meanArray);
        /**
        * Helper routine to center an array (mean for each row=0).
        * @param[in] inputData row major data to update.
        */
        static void centerRow2DData(dblarray &inputData);
        /**
        * Helper routine to center an array (mean for each row=0).
        * @param[in] inputData row major data to update.
        * @param[in] MeanName string with filename to save fits on.
        */
        static void centerRow2DData(dblarray &inputData,char *MeanName);
        /**
        * Helper routine to project an array on l2 ball (row l2norm <1).
        * @param[in] inputData row major data to update.
        */
        static void l2ProjRow2DData(dblarray &inputData);
        /**
        * Helper routine to center and project an array on l2 ball
            (row l2norm <1).
        * @param[in] inputData row major data to update.
        */
        static void l2NormalizeRow2DData(dblarray &inputData);
        /**
        * Helper routine to compute the Gram Matrix [Ncols x Ncols].
        * @param[in] mat matrix row to compute the Gram matrix.
        */
        static dblarray* gram(dblarray &mat);
        /**
         * Helper routine to compute the Gram Matrix after normalizing columns
           in l2 norm (except 0 column) [Ncols x Ncols].
          * @param[in] mat matrix row to compute the Gram matrix.
        */
        static dblarray* unitL2Gram(dblarray &mat);
        /**
         * Helper routine to get spectral radius of a matrix.
          * @param[in] Matrix matrix to get spectral radius.
          * @param[in] xsubi random seed for initialization.
          * @param[in] nit_max maximal number of iterations.
          * @return the spectral radius
        */
        template<typename PARAM>
        static double getMatrixL2Norm(to_array<PARAM,true> * Matrix,
                                        unsigned short xsubi[3], int nit_max);

        /**
        * Helper routine to multiply two row major matrices.
        * @param[in] mat1 matrix1 [Ncol1,Nrow1].
        * @param[in] mat2 matrix2 [Ncol2,Nrow2].
        * @return the product of mat1 and mat2 [Ncol2,Nrow1] (Ncol1=Nrow2)
        */
        static dblarray* matMult(dblarray &mat1, dblarray &mat2);
        /**
        * Helper routine to multiply two row major matrices.
        * @param[in] Mout matrix to contain the output [Ncol2,Nrow1].
        * @param[in] M1 matrix1 [Ncol1,Nrow1].
        * @param[in] M2 matrix2 [Ncol2,Nrow2].
        * @return the product of mat1 and mat2 [Ncol2,Nrow1] (Ncol1=Nrow2)
        *@note if Mout is defined (not a null pointer), it uses the same object.
          Otherwise the object is created
        */
        static void matMult(dblarray* &Mout, dblarray &mat1, dblarray &mat2);
        static void matMult(dblarray &Mout, dblarray &mat1, dblarray &mat2);

        /**
        * Helper routine to multiply a matrix by the transpose of another.
        * @param[in] mat1 matrix1 [Ncol1,Nrow1].
        * @param[in] mat2 matrix2 [Ncol2,Nrow2].
        * @return the product of mat1 and transpose(mat2) [Nrow2,Nrow1]
        */
        static dblarray* matMultTransp(dblarray &mat1, dblarray &mat2);
        /**
        * Helper routine to multiply a matrix by the transpose of another.
        * @param[in] Mout matrix to contain the output [Ncol1,Nrow1].
        * @param[in] M1 matrix1 [Ncol1,Nrow1].
        * @param[in] M2 matrix2 [Ncol2,Nrow2].
        * @return the product of mat1 and transpose(mat2) in Mout [Nrow2,Nrow1]
        *@note if Mout is defined (not a null pointer), it uses the same object.
            Otherwise the object is created
        * @brief Multiply M1 * tranpose(M2), both row major:
                Mout[kj+ki*Ny1]=sum_kk M1[kk+ki*Nx1]*M2[kk+kj*Nx2]
        */
        static void matMultTransp(dblarray* &Mout, dblarray &M1, dblarray &M2);
        static void matMultTransp(dblarray &Mout, dblarray &mat1, dblarray &mat2);
        /**
        * Helper routine to multiply a matrix by another and then transpose.
        * @param[in] mat1 matrix1 [Ncol1,Nrow1].
        * @param[in] mat2 matrix2 [Ncol2,Nrow2].
        * @return the tranposed product of mat1 and mat2 [Nrow1,Ncol2]
        */
        static dblarray* transpMatMult(dblarray &mat1, dblarray &mat2);
        /**
        * Helper routine to perform matrix vector multiplication.
        * @param[in] mat matrix [Ncol1,Nrow1].
        * @param[in] vect matrix of vectors [Ncol1,Nrow2].
        * @param[in] index index of the vector for multiplication (>0 <Nrow2).
        * @return the matrix product [Nrow1]
        */
        static dblarray* matRowVec(const dblarray &mat, const dblarray &vect,
                                                            unsigned int index=0);
        /**
        * Helper routine to perform inner product.
        * @param[in] vector1 first vector .
        * @param[in] vector2 second vector.
        * @param[in] stride1 stride for first vector.
        * @param[in] stride2 stride for second vector.
        * @param[in] Nelems vector dimension.
        * @param[in] metric metric to use in the inner product.
        * @param[in] useMetric use a metric [NelemsxNelems] or [Nelems].
        * @param[in] diagMetric use a diagonal metric [size Nelems].
        * @return the inner product
        */
        static double dotProduct(const double* vector1, const double* vector2,
                            size_t stride1,size_t stride2,size_t Nelems,
                            const double* metric,bool useMetric,bool diagMetric);
        /**
        * Helper routine to perform inner product.
        * @param[in] vector1 first vector .
        * @param[in] vector2 second vector.
        * @param[in] metric metric to use in the inner product.
        * @param[in] useMetric use a metric [NelemsxNelems] or [Nelems].
        * @param[in] diagMetric use a diagonal metric [size Nelems].
        * @return the inner product
        */
        static double dotProduct(const dblarray  &vector1, const dblarray &vector2,
                            const dblarray &metric,bool useMetric,bool diagMetric);
        //pure virtual processing functions to implement
        /**
        * Interface routine to perform cholesky decomposition product.
        * @param[in] PosDefMat positive definite matrix.
        * @return the cholesky decomposition
        */
        virtual dblarray* cholesky(dblarray& PosDefMat)=0;
        /**
        * Interface routine to perform svd.
        * @param[in] Mat matrix to decompose.
        * @param[out] U left unitary matrix of singular vectors.
        * @param[out] S diagonal with singular values.
        * @param[out] Vt right unitary matrix of singular vector.
        * @return the singular value decomposition, Mat=U S Vt
        */
        virtual void svd(dblarray& Mat, dblarray& U, dblarray& S,
                                                                dblarray& Vt)=0;
        /**
        * Interface routine to perform inverse of a square matrix.
        * @param[in] SquareMat square full rank matrix to invert.
        * @return the inverse
        */
        virtual dblarray* invert(dblarray& SquareMat)=0;
        /**
        * Interface routine to perform pseudo inverse.
        * @param[in] Mat input matrix.
        * @param[in] thr relative threshold for eigenvalues [truncation].
        * @return the pseudo inverse of the matrix
        */
        virtual dblarray* pinvert(dblarray& Mat, double thr)=0;
        /**
        * Interface routine to perform pseudo inverse.
        * @param[in] Mat input matrix.
        * @return the pseudo inverse of the matrix
        * @brief it uses a max inverse condition number
          @see MatrixDecomposition::m_invcondnb
        */
        virtual dblarray* pinvert(dblarray& Mat){
                                             return pinvert(Mat, m_invcondnb);}

        dblarray* sqrtSVD(dblarray& Mat, double relthr=0);

      virtual void solve_UT(dblarray& A, dblarray& b, dblarray& x);
      virtual void solve_LT(dblarray& A, dblarray& b, dblarray& x);
   private:

};

/**
 * @class SapDecomposition
 * @brief Class implementing decompositions using the default Sap package.
*/
class SapDecomposition: public MatrixDecomposition {
    protected:
    public:
        //constructors/destructors
        /**
        * Standart constructor.
        */
        SapDecomposition(){}
        /**
        * Standart destructor.
        */
        ~SapDecomposition(){}
        //setters
        //getters
        /**
        * Get the signature of the object.
        * @return signature of SapDecomposition.
        */
        virtual std::string getSignature() const {return std::string("Sapdev");}
        //processing function
        /**
        * Implementation of cholesky decomposition product.
        * @param[in] PosDefMat positive definite matrix.
        * @return the cholesky decomposition, NOT IMPLEMENTED YET
        */
        virtual dblarray* cholesky(dblarray& PosDefMat);
        /**
        * Implementation of inverse of a square matrix.
        * @param[in] SquareMat square full rank matrix to invert.
        * @return the inverse
        */
        virtual dblarray* invert(dblarray& SquareMat);
        /**
        * Implementation of svd.
        * @param[in] Mat matrix to decompose.
        * @param[out] U left unitary matrix of singular vectors.
        * @param[out] S diagonal with singular values.
        * @param[out] Vt right unitary matrix of singular vector.
        * @return the singular value decomposition, Mat=U S Vt
        */
        virtual void svd(dblarray& Mat, dblarray& U, dblarray& S,
                                                                dblarray& Vt);
        /**
        * Implementation of pseudo inverse.
        * @param[in] Mat input matrix.
        * @param[in] relthr relative threshold for eigenvalues [truncation].
        * @return the pseudo inverse of the matrix
        */
        virtual dblarray* pinvert(dblarray& Mat,double relthr);
};

#ifdef _ATLAS_DEC
/**
 * @class AtlasDecomposition
 * @brief Class implementing decompositions using ATLAS,LAPACK(+E) libraries.
*/
class AtlasDecomposition: public MatrixDecomposition {
    public:
        //constructors/destructors
        /**
        * Standart constructor.
        */
        AtlasDecomposition(){}
        /**
        * Standart destructor.
        */
        ~AtlasDecomposition(){}
        //getters
        /**
        * Get the signature of the object.
        * @return signature of AtlasDecomposition.
        */
        virtual std::string getSignature() const {return std::string("Atlas");}
        //processing function
        /**
        * Implementation of cholesky decomposition product.
        * @param[in] PosDefMat positive definite matrix.
        * @return the cholesky decomposition, NOT IMPLEMENTED YET
        */
        virtual dblarray* cholesky(dblarray& PosDefMat);
        /**
        * Implementation of inverse of a square matrix.
        * @param[in] SquareMat square full rank matrix to invert.
        * @return the inverse
        */
        virtual dblarray* invert(dblarray& SquareMat);
        /**
        * Implementation of svd.
        * @param[in] Mat matrix to decompose.
        * @param[out] U left unitary matrix of singular vectors.
        * @param[out] S diagonal with singular values.
        * @param[out] Vt right unitary matrix of singular vector.
        * @return the singular value decomposition, Mat=U S Vt
        */
        virtual void svd(dblarray& Mat, dblarray& U, dblarray& S,
                                                                dblarray& Vt);
        /**
        * Implementation of pseudo inverse.
        * @param[in] Mat input matrix.
        * @param[in] relthr relative threshold for eigenvalues [truncation].
        * @return the pseudo inverse of the matrix
        */
        virtual dblarray* pinvert(dblarray& Mat,double relthr);
};
#endif

#ifdef _ARMA_DEC
/**
 * @class ArmadilloDecomposition
 * @brief Class implementing decompositions using ARMADILLO,LAPACK(+E) libraries.
*/
class ArmadilloDecomposition: public MatrixDecomposition {

    protected:
        ArmadilloSVD  m_svdMethod;

    public:
        //constructors/destructors
        /**
        * Standart constructor.
        */
        ArmadilloDecomposition(){m_svdMethod=DC;}
        /**
        * Constructor specifying the type of method for svd.
        * @param[in] svdMeth method to compute SVD ("d"c, divide and conquer,
            vs "std" standard method)
        */
        ArmadilloDecomposition(ArmadilloSVD svdMeth){m_svdMethod= svdMeth;}
        /**
        * Standart destructor.
        */
        ~ArmadilloDecomposition(){}
        //setters
        /**
        * Specify the type of method for svd.
        * @param[in] svdMeth method to compute SVD ("d"c, divide and conquer,
            vs "std" standard method)
        */
        void setSvdMethod(const ArmadilloSVD svdmethod){m_svdMethod=svdmethod;}
        //getters
        /**
        * Get the method to compute SVD.
        * @return signature of method for SVD computation.
        */
        ArmadilloSVD getSvdMethod() const {return m_svdMethod;}
        /**
        * Get the signature of the object.
        * @return signature of ArmadilloDecomposition.
        */
        virtual std::string getSignature() const {
            if(m_svdMethod==DC) return std::string("Armadillo");
            else return std::string("ArmadilloSTD");
        }
        //processing functions
        /**
        * Implementation of cholesky decomposition product.
        * @param[in] PosDefMat positive definite matrix.
        * @return the cholesky decomposition, NOT IMPLEMENTED YET
        */
        virtual dblarray* cholesky(dblarray& PosDefMat);
        /**
        * Implementation of inverse of a square matrix.
        * @param[in] SquareMat square full rank matrix to invert.
        * @return the inverse
        */
        virtual dblarray* invert(dblarray& SquareMat);
        /**
        * Implementation of svd.
        * @param[in] Mat matrix to decompose.
        * @param[out] U left unitary matrix of singular vectors.
        * @param[out] S diagonal with singular values.
        * @param[out] Vt right unitary matrix of singular vector.
        * @return the singular value decomposition, Mat=U S Vt
        */
        virtual void svd(dblarray& Mat, dblarray& U, dblarray& S,dblarray& Vt);
        /**
        * Implementation of pseudo inverse.
        * @param[in] Mat input matrix.
        * @param[in] relthr relative threshold for eigenvalues [truncation].
        * @return the pseudo inverse of the matrix
        */
        virtual dblarray* pinvert(dblarray& Mat,double thr);
        //BEWARE FOR THIS LAST FUNCTION, THIS IS AN ABSOLUTE THRESHOLD (NOT
        //RELATIVE AS FOR OTHER DECOMPOSITIONS)
        /**
        * Implementation of solve Ax=b, with A Upper/Lower-Triangular
        */
        virtual void solve_UT(dblarray& A, dblarray& b, dblarray& x);
        virtual void solve_LT(dblarray& A, dblarray& b, dblarray& x);
};
#endif

#include "MatrixL2Norm.tcc"

#endif //MatrixDecomposition_H

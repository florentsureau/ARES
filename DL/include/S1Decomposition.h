/******************************************************************************
 **                   Copyright (C) 2017 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Florent Sureau
 **
 **    Date:  01-02/2017
 **
 **    File:  S1Decomposition.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces classes for calculus on (S1)^n.
 **    -----------
 ******************************************************************************/
 

#ifndef S1Decomposition_H
#define S1Decomposition_H
#include "MatrixDecomposition.h"
#include "DLException.h"
#define EpsAngle 1e-10 //"0" angle difference
//Note: dblarray is row major ordering of the buffer if(x,y)=(col,row).
/**
 * @class S1Decomposition 
 * @brief Class essentially implementing the basic helper routines for S1 data
          parameterized by angle.
*/
class S1Decomposition
{
    protected:
        size_t m_Nx;/**< Number of S1 measurements (S1.*/
        size_t m_Ny;/**< Number of rows.*/
        int m_naxis;/**< Number of dimensions.*/
    public:
        /**
         * Standart constructor.
        */
        S1Decomposition():m_Nx(0), m_Ny(0),m_naxis(0){}
        /**
         * Standart destructor.
         * @return number of DL iterations.
        */
        virtual ~ S1Decomposition(){}
        //members
        //setters
        //processing helpers
        /**
         * Wrap real angle in [-pi,pi[
          * @param[in] input Angle to wrap
          * @return wrapped angle in [-pi,pi[
        */
        static inline double wrap_real_2pi(double input){
            int N2PI=floor((floor(input/M_PI)+1)/2);
            return input-N2PI*2*M_PI;
        }

        /**
         * Compute S1 Distance between two points in (S1)^Nels.
          * @param[in] S1vec1 first point in (S1)^Nels. 
          * @param[in] S1vec2 second point in (S1)^Nels
          * @param[in] Nels number of coordinates in vectors
          * @return distance: sum of squared wrapped distance for coordinates 
          * @brief: inputs are assumed contiguous
        */
        static double computeS1Dist(double* S1vec1,double* S1vec2, size_t Nels);
        /**
         * Compute S1 Distance between two points in (S1)^Nels.
          * @param[in] S1vec1 first point in (S1)^Nels. 
          * @param[in] stride1 stride for first component. 
          * @param[in] S1vec2 second point in (S1)^Nels
          * @param[in] stride2 stride for second component. 
          * @param[in] Nels number of coordinates in vectors
          * @return distance: sum of squared wrapped distance for coordinates 
        */
        static double computeS1Dist(double* S1vec1,size_t stride1,double* S1vec2,
                                                   size_t stride2, size_t Nels);

        /**
         * Compute linfNorm on logmap between two points in (S1)^Nels .
          * @param[in] S1vec1 first point in (S1)^Nels. 
          * @param[in] S1vec2 second point in (S1)^Nels
          * @param[in] Nels number of coordinates in vectors
          * @return distance: max of wrapped l1 norm on logmap for coordinates 
          * @brief: inputs are assumed contiguous
        */
        static double computeS1linf(double* S1vec1,double* S1vec2,size_t Nels);
        /**
         * Compute linfNorm on logmap between two points in (S1)^Nels .
          * @param[in] S1vec1 first point in (S1)^Nels. 
          * @param[in] stride1 stride for first component. 
          * @param[in] S1vec2 second point in (S1)^Nels
          * @param[in] stride2 stride for second component. 
          * @param[in] Nels number of coordinates in vectors
          * @return distance: max of wrapped l1 norm on logmap for coordinates 
        */
        static double computeS1linf(double* S1vec1,size_t stride1,double* S1vec2,
                                                   size_t stride2, size_t Nels);

        /**
         * Check second point in the max strongly convex ball centered 
            on first point .
          * @param[in] S1vec1 first point in (S1)^Nels. 
          * @param[in] S1vec2 second point in (S1)^Nels
          * @param[in] Nels number of coordinates in vectors
          * @return true if the logmap per coordinate is in the S1 ball of radius 
            PI/2
          * @brief: inputs are assumed contiguous
        */
        static bool checkMaxCvxBall(double* S1vec1,double* S1vec2, size_t Nels);
        /**
         * Check second point in the max strongly convex ball centered 
            on first point .
          * @param[in] S1vec1 first point in (S1)^Nels. 
          * @param[in] stride1 stride for first component. 
          * @param[in] S1vec2 second point in (S1)^Nels
          * @param[in] stride2 stride for second component. 
          * @param[in] Nels number of coordinates in vectors
          * @return true if the logmap per coordinate is in the S1 ball of radius 
            PI/2
        */
        static bool checkMaxCvxBall(double* S1vec1,size_t stride1,double* S1vec2,
                                                   size_t stride2, size_t Nels);

        /**
         * Add an angle per row of input array.
          * @param[in/out] input input angle array in [-pi,pi[, where in place 
            addition will be performed
          * @param[in] RowAdd angle array with one element per row of input array
        */
        static void addAngleRow(dblarray* input, const dblarray& RowAdd);

        /**
         * Compute Frechet function assuming input sorted.
          * @param[in] testPoint reference. 
          * @param[in] Elements points to compute for Frechet function, sorted
            in [-pi, pi[
          * @param[in] Nels number of points
          * @return Frechet function: sum of squared distances from the points 
            to the reference point
        */
       static double computeSortedFrechet(double testPoint,
                                                double* Elements, size_t Nels);
        /**
         * Compute Frechet function first sorting the input.
          * @param[in] testPoint reference. 
          * @param[in] Elements points to compute for Frechet function
          * @param[in] Nels number of points
          * @return Frechet function: sum of squared distances from the points 
            to the reference point
        */
       static double computeUnsortedFrechet(double testPoint,
                                                double* Elements, size_t Nels);
        /**
         * Compute Frechet function assuming input sorted.
          * @param[in] testPoint reference. 
          * @param[in] Elements points to compute for Frechet function, sorted
            in [-pi, pi[
          * @param[in] Weights weights (theoretically summing to 1) of each coef
          * @param[in] Nels number of points
          * @return Frechet function: sum of weighted squared distances from the points 
            to the reference point
        */
       static double computeSortedWeightedFrechet(double testPoint,
                                double* Elements, double* Weights, size_t Nels);
        /**
         * Compute Frechet function Grad Norm assuming input sorted.
          * @param[in] testPoint reference. 
          * @param[in] Elements points to compute for Frechet function, sorted
            in [-pi, pi[
          * @param[in] Weights weights of each element
          * @param[in] Nels number of points
          * @param[in] rightLim compute right limit instead of left limit (i.e
            consider angles in ]-pi,pi].
          * @return Frechet function: sum of squared weighted difference from the points 
            to the reference point
        */
       static double computeSortedWeightedFrechetGradNorm(double testPoint,
            double* Elements, double* Weights, size_t Nels,bool rightLim=false);


        /**
         * Compute Frechet Mean.
          * @param[in] S1Angle buffer where to get input data in [-pi,pi[. 
          * @param[out] Mean buffer where to store Frechet Mean
          * @param[in] Nx number of S1 elements to average for each mean
          * @param[in] Nvec number of average to compute
          
        */
        static void computeFrechetMean(double* S1Angle,double* Mean,
                                                        size_t Nx,size_t Nvec);
        /**
         * Compute weighted Frechet Mean.
          * @param[in] S1Angle buffer where to get input data in [-pi,pi[. 
          * @param[in] weights weights associated to each angle. 
          * @param[out] Mean buffer where to store weighted Frechet Mean
          * @param[in] Nx number of S1 elements to average for each mean
          * @param[in] Nvec number of average to compute
          
        */
        static void computeWeightFrechetMean(double* S1Angle,double* weights,
                        double* WFrechMean,size_t Nx,size_t Nvec);
        /**
         * Compute min of l2 norm of weighted Frechet gradient.
          * @param[in] S1Angle buffer where to get input data in [-pi,pi[. 
          * @param[in] weights weights associated to each angle. 
          * @param[out] WFrechGradMin buffer where to store min
          * @param[in] Nx number of S1 elements to average for each mean
          * @param[in] Nvec number of average to compute
          
        */
        static void computeWeightFrechetMinGradNorm(double* S1Angle,
                  double* weights, double* WFrechGradMin,size_t Nx,size_t Nvec);

        /**
         * Compute Frechet Mean.
          * @param[in] S1Angle array where we compute the mean per row.
          * @return array of Frechet mean for each row.
        */
        static dblarray* computeFrechetMean(const dblarray& S1Angle);


        /**
         * Compute Weighted Frechet Mean.
          * @param[in] S1Angle array where we compute the mean per row.
          * @return array of Frechet mean for each row.
        */
        static dblarray* computeWeightedFrechetMean(const dblarray& S1Angle, 
            const dblarray& Weights);
            
            

        /**
         * subtract Frechet Mean.
          * @param[in] S1Angle array where we subtract to each row its mean.
          * @return array of Frechet mean for each row.
        */
        static dblarray* subtractFrechetMean(dblarray& S1Angle);

        /**
         * subtract Frechet Mean for column data.
          * @param[in] S1Angle array where we subtract to each column its mean.
          * @return array of Frechet mean for each column.
        */
        static dblarray* subtractFrechetMeanCol(dblarray& S1Angle);

        /**
         * Compute weighted inputs 2pi shifted according to reference:
           Compute W_{ik}=sum_j weight{ij} (input_{ik}+2pi*S_{ijk})
           , where (S_{ijk})in {0,-1,1} and input_{ik}+2pi*S_{ijk})-Ref_{jk} in
           [-pi,pi[
          * @param[in] Input array of angle values (pixels as columns, 
            samples as rows).
          * @param[in] Reference array of angle values (references as columns, 
            pixels as rows).
          * @param[in] weighting array (referencs as ,columns samples as rows).
          * @param[in] shiftWInput ptr to an array  where the shifted 
                input will be stored.
          * @return true if input need to be shifted again.
        */
        static bool weightShiftInputToReference(const dblarray& Input,
                const dblarray& Reference,const dblarray& weighting, 
                                                        dblarray* &shiftWInput);
                                                        
        /**
         * Compute Affine weighted residual:
           Compute W_{ik}=sum_j weight{ij} (input_{ik}+2pi*S_{ijk}-Atoms_{jk})
           , where (S_{ijk})in {0,-1,1} and input_{ik}+2pi*S_{ijk})-Atoms_{jk} in
           [-pi,pi[
          * @param[in] Input array of angle values (pixels as columns, 
            samples as rows).
          * @param[in] Atoms array with atoms as columns, pixels as rows.
          * @param[in] weighting array (referencs as ,columns samples as rows).
          * @param[in] Residual ptr to an array  where the residual is stored 
        */
        void computeAffWResidual(const dblarray& Input,const dblarray& Atoms,
                                   const dblarray& Weight, dblarray* & WResidual);

        /**
         * Compute Log Euclidean residual:
           Compute W_{ik}=sum_j weight{ij}(2pi*S_{ijk}+Atoms_{jk} -Tgt_{jk})
             -input_{ik}), where (S_{ijk})in {0,-1,1} and 2pi*S_{ijk}+Atoms_{jk}
            -Tgt_{ik} (or -Tgt_{i}) in [-pi,pi[
          * @param[in] Input array of angle values (pixels as columns, 
            samples as rows).
          * @param[in] TgtPlane array of angle values used as reference for 
              tangent plane(either of size number of samples, or 
              number of pixelsxnumber of samples).
          * @param[in] Atoms array with atoms as columns, pixels as rows.
          * @param[in] weighting array (referencs as ,columns samples as rows).
          * @param[in] Residual ptr to an array  where the residual is stored 
        */
        void computeLogEuclResidual(const dblarray& Input,
                                const dblarray& TgtPlane, const dblarray& Atoms,
                                const dblarray& Weight, dblarray* & WResidual);

        /**
         * Compute gradient of weighted l2 norm on residual:
           Compute W_{jk}=sum_i w_{ij} sum_j' weight{ij'} (input_{ik}+
                                                       2pi*S_{ij'k}-Atoms_{j'k})
           , where (S_{ij'k})in {0,-1,1} and input_{ik}+2pi*S_{ij'k})-Ref_{j'k} 
            is in [-pi,pi[
          * @param[in] Input array of angle values (pixels as columns, 
            samples as rows).
          * @param[in] Atoms array with atoms as columns, pixels as rows.
          * @param[in] weighting array (referencs as ,columns samples as rows).
          * @param[in] GradL2 ptr to an array  where the residual is stored 
        */
        void computeGradientAffWL2Norm(const dblarray& Input,
               const dblarray& Atoms,const dblarray& Weight, dblarray* &GradL2);

        /**
         * Compute gradient of weighted l2 norm on residual:
           Compute W_{jk}=sum_i w_{ij} (sum_j' weight{ij'}  2pi*S_{ij'k}+
                Atoms_{j'k}-Tgt_{ik}-input_{ik})
           , where (S_{ij'k})in {0,-1,1} and 2pi*S_{ij'k}+Atoms_{j'k}- Tgt_{ik}
            (or -Tgt_{i}  is in [-pi,pi[
          * @param[in] Input array of angle values (pixels as columns, 
            samples as rows).
          * @param[in] TgtPlane array of angle values used as reference for 
              tangent plane(either of size number of samples, or 
              number of pixelsxnumber of samples).
          * @param[in] Atoms array with atoms as columns, pixels as rows.
          * @param[in] weighting array (referencs as ,columns samples as rows).
          * @param[in] GradL2 ptr to an array  where the residual is stored 
        */
        void computeGradientLogEuclL2Norm(const dblarray& Input,
                                const dblarray& TgtPlane,const dblarray& Atoms,
                                    const dblarray& Weight, dblarray* &GradL2);

        /**
         * Compute gradient of weighted l2 norm on residual:
           Compute W_{jk}=sum_i w_{ij} sum_j' weight{ij'} (input_{ik}+
                                                       2pi*S_{ij'k}-Atoms_{j'k})
           , where (S_{ij'k})in {0,-1,1} and input_{ik}+2pi*S_{ij'k})-Ref_{j'k} 
            is in [-pi,pi[
          * @param[in] WResidual ptr to an array  where are weighted residuals
          * @param[in]  Weight array (referencs as ,columns samples as rows).
          * @param[in] GradL2 ptr to an array  where the residual is stored 
        */
        void computeGradientAffWL2Norm(const dblarray& WResidual, 
                                    const dblarray& Weight, dblarray* &GradL2);


        /**
         * Compute gradient of weighted l2 norm on residual:
           Compute W_{jk}=sum_i w_{ij} (input_{ik}-sum_j' weight{ij'} 
                               2pi*S_{ij'k}+Atoms_{j'k}-Tgt_{ik} (or -Tgt_{i})
           , where (S_{ij'k})in {0,-1,1} and 2pi*S_{ij'k}+Atoms_{j'k}- Tgt_{ik}
            (or -Tgt_{i}  is in [-pi,pi[
          * @param[in] WResidual ptr to an array  where are weighted residuals
          * @param[in]  Weight array (referencs as ,columns samples as rows).
          * @param[in] GradL2 ptr to an array  where the residual is stored 
        */
        void computeGradientLogEuclL2Norm(const dblarray& WResidual, 
                                    const dblarray& Weight, dblarray* &GradL2);

   private:
   
};


#endif //S1Decomposition_H

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
 **    File:  MatrixL2Norm.tcc
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Implementation of Templatized Power Method.
 **    -----------
 ******************************************************************************/


 template<typename PARAM>
 double MatrixDecomposition::getMatrixL2Norm(to_array<PARAM,true>* Matrix,
                                          unsigned short xsubi[3],int nit_max) {
     FINFO_;
     //Power Method
     //see Estimating the matrix p-norm, NJ Higham, Numer. Math. 62:539-555, 1992
     //Nx=Natoms, Ny=Npixels
     long Nx=Matrix->nx(), Ny=Matrix->ny();
     double init_vector[Nx],vec_fwd[Ny],vec_bwd[Nx],vec_bwd2[Nx];
     long kl,kl1,kl2, offset;
     double gamma=-1.;
     double Nrm2In=0.,Nrm2fwd,Nrm2bwd,Totalbwd;
     PARAM *buffer_mat= Matrix->buffer();

     for(kl=0;kl<Nx;kl++) {
         init_vector[kl]=erand48(xsubi);
         Nrm2In+=init_vector[kl]*init_vector[kl];
     }

     Nrm2In =sqrt(Nrm2In);
     for(kl=0;kl<Nx;kl++) init_vector[kl] /= Nrm2In;

     for(int kit=0;kit < nit_max;kit++) {
         Nrm2fwd=0.;
         #ifdef _OPENMP
             #pragma omp parallel for default(none) reduction(+: Nrm2fwd)\
                             shared(vec_fwd,init_vector, buffer_mat,Nx,Ny)\
                             private(kl1,kl2,offset) schedule(static)
         #endif
         for(kl1=0;kl1<Ny;kl1++) {
             vec_fwd[kl1]=0.;
             offset=kl1* Nx;
             for(kl2=0;kl2<Nx;kl2++) {
                 vec_fwd[kl1] += (double) buffer_mat[offset+kl2] * init_vector[kl2];
             }
             Nrm2fwd+= vec_fwd[kl1]* vec_fwd[kl1];
        }
         Nrm2fwd =sqrt(Nrm2fwd);
         for(kl1=0;kl1<Ny;kl1++) vec_fwd[kl1] /= Nrm2fwd;

         for(kl1=0;kl1<Nx;kl1++) vec_bwd[kl1]=0.;
         #ifdef _OPENMP
             #pragma omp parallel for default(none) \
                                     shared(vec_bwd, vec_fwd, buffer_mat,Nx,Ny)\
                                             private(kl1,kl2) schedule(static)
         #endif
         for(kl1=0;kl1<Nx;kl1++) {
             for(kl2=0;kl2<Ny;kl2++)
                     vec_bwd[kl1]+=(double) buffer_mat[kl2* Nx+kl1]* vec_fwd[kl2];
         }
         Nrm2bwd =0.;
         for(kl1=0;kl1<Nx;kl1++) Nrm2bwd += vec_bwd[kl1]* vec_bwd[kl1];
         Nrm2bwd =sqrt(Nrm2bwd);
         Totalbwd=0.;
         for(kl1=0;kl1<Nx;kl1++) {
             vec_bwd2[kl1] = vec_bwd[kl1]* init_vector[kl1];
             Totalbwd+=vec_bwd2[kl1];
         }
         if(Nrm2bwd < Totalbwd) {
             gamma= Nrm2fwd;
             break;
         }
         for(kl1=0;kl1<Nx;kl1++) init_vector[kl1] = vec_bwd[kl1] / Nrm2bwd;
     }
     /*if(gamma < 0.)
         printf("Spectral Radius not converged in %d iterations\n",nit_max);
     printf("gamma =%lf\n", Nrm2fwd);*/
     gamma = Nrm2fwd;
     FOFNI_;
     return gamma;
 }

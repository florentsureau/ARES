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
 **    File:  gslRNG.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces classes for Random number generation using GSL.
 **    -----------
 **
 **
 ******************************************************************************/


#ifndef GSLRNG_H
#define GSLRNG_H
#include "TempArray.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <time.h>

class gslRNG {
   public:
      //Constructors/Destructors
      gslRNG(): m_rng(NULL), m_perm(NULL){initRNG();}
      gslRNG(const gsl_rng_type* T): m_rng(NULL), m_perm(NULL){initRNG(T);}
      virtual ~gslRNG() {
         if(m_rng!=NULL) gsl_rng_free(m_rng);
         if(m_perm!=NULL) gsl_permutation_free(m_perm);
      }
      void initRNG(const gsl_rng_type* T=gsl_rng_default){
         if(m_rng!=NULL) gsl_rng_free (m_rng);
         m_rng =gsl_rng_alloc(T);
      }
      void createSeed(){
         // Generating a random seed using  /dev/urandom
         m_seed=(size_t) time(NULL);
         setSeed(m_seed);
      }
      //Setters
      void setSeed(size_t seed){
         m_seed=seed;
         gsl_rng_set(m_rng, m_seed);
      }
      //Getters
      size_t getSeed() const {return m_seed;}
      std::string getSignature(){return std::string("gslRNG");}

      //Processing Functions
      gsl_permutation* createPermutation(size_t Nelem){
         if(m_rng ==NULL) initRNG();
         if(m_perm!=NULL) gsl_permutation_free(m_perm);
         m_perm =gsl_permutation_calloc (Nelem);
         if(m_perm ==NULL){
            std::cout<<"Cannot Allocate Memory for permutation:"
                                                            <<Nelem<<std::endl;
            exit(EXIT_FAILURE);
         }
         gsl_ran_shuffle (m_rng, m_perm->data, Nelem, sizeof(size_t));
         return m_perm;
      }
      size_t getUniformInt(size_t maxVal){
         if(m_rng ==NULL) initRNG();
         return gsl_rng_uniform_int(m_rng,maxVal);
      }
      double getUniformDbl(double maxVal){
         if(m_rng ==NULL) initRNG();
         return gsl_rng_uniform(m_rng)*fabs(maxVal);
      }
      size_t getBernoulli(double psuccess){
         if(m_rng ==NULL) initRNG();
         return gsl_ran_bernoulli(m_rng,psuccess);
      }
      double getNormal(double mean, double sigma){
         if(m_rng ==NULL) initRNG();
         return (gsl_ran_gaussian(m_rng,sigma)+mean);
      }

   protected:
      gsl_rng* m_rng;
      gsl_permutation* m_perm;
      size_t m_seed;

   private:
};
#endif

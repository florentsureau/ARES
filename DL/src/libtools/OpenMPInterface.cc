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
 **    Date:  09/2017
 **
 **    File:  OpenMPInterface.cpp
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Interface to set OpenMP threads
 **    -----------
 ** - 09/2017: Added automatic/manual setting
 ******************************************************************************/

#include "OpenMPInterface.h"

unsigned int OMPThreadHandling::mInnerThreads=1;
unsigned int OMPThreadHandling::mOuterThreads=1;
unsigned int OMPThreadHandling::mTotalThreads=1;

using namespace std;

void OMPThreadHandling::getStatus() {

   #ifdef _OPENMP
	    cout << "Application was compiled with OpenMP support," << endl;
    	if (omp_get_max_threads() == 1) cout << "but running with one process only." << endl;
    	else  cout << "running with up to " << omp_get_max_threads()
           					<< " processes." << endl;
        if (omp_get_nested()){
           	cout << "Nested Parallelism allowed" <<endl;
           	cout << "running outer loops with " << mOuterThreads
           					<< " processes." << endl;
           	cout << "running inner loops with " << mInnerThreads
           					<< " processes." << endl;

         } else cout << "running with" << mOuterThreads
           					<< " processes." << endl;
    #else
    	 cout << "Application was compiled without OpenMP support;" << endl
        	 << "running in scalar mode." << endl;
	mInnerThreads=mOuterThreads=1;
    #endif
  };


void OMPThreadHandling::setNThreads(const int outer, const int inner) {
     #ifdef _OPENMP
 	int nthr=omp_get_max_threads();

	if(outer * inner > nthr) {
		cout << "Too many threads required (" <<outer << "*" << inner <<" vs "<< nthr <<" procs)"<<endl;
		exit(EXIT_FAILURE);
	} else {
		if(outer <= 0) mOuterThreads=1; else mOuterThreads=(unsigned int) outer;
		if(inner <=0)  mInnerThreads=1; else mInnerThreads=(unsigned int) inner;
	}
	mTotalThreads=mOuterThreads*mInnerThreads;
	omp_set_num_threads(mTotalThreads);
        std::cout<<"SET THREADS"<<omp_set_num_threads<<std::endl;
    #endif

  };

void OMPThreadHandling::autoSetThreads(const int outer_max,const int inner_max,const int nb_procs) {
	//Simple rules for setting the number of outer/inner threads:
	//First set the maximal number of procs available if not set explicitely:
	//1) if nprocs  < 4, then we use all of them
	//2) else use only half of them (we do not need so much computing power)
	//Then dispatch the procs in outer loops and inner loops knowing that:
	//outer_max (inner_max) is the maximal number of threads useful for outer (inner) loops, set to the nb of procs available  if equal to 0.
	//The criterion is the following:
	//Maximize the number of procs used, and if several combinations, maximize the number of proc for outer loops threads (where we expect more computing time gain)
  	int nbprocs_avail,nthr,outer_l=outer_max,inner_l=inner_max;
  	omp_set_nested(1);

	nbprocs_avail=omp_get_max_threads();

	if(nb_procs==0) {
		if(nbprocs_avail<4) omp_set_num_threads(nbprocs_avail);
		else omp_set_num_threads((int) (nbprocs_avail/2)); //Use only 50% of procs
	} else {
		if(nb_procs > nbprocs_avail) 	{
			cout << "Too many procs required ("<<nbprocs_avail<<" required vs "<<nb_procs <<" avail)" <<endl ;
			exit(EXIT_FAILURE);
		} else {
			nbprocs_avail=nb_procs;
			omp_set_num_threads((int) (nbprocs_avail));
		}
	}
	nthr=omp_get_max_threads();
	if(!omp_get_nested()) setNThreads(nthr,1);
	else {
		if((outer_max > nthr) || (outer_max<=1)) outer_l=nthr;
		if((inner_max > nthr) || (inner_max<=1)) inner_l=nthr;

		int k,divisor,outer,inner;
		int MaxPotRed=0;
		for(k=outer_l-1;k>=0;k--) {
			divisor=(int) floor(nthr /(k+1));
			if((divisor <= inner_l) && (divisor*(k+1) > MaxPotRed)) {
				MaxPotRed=divisor*(k+1);
				outer=k+1;
				inner=divisor;
			}
		}
		setNThreads(outer,inner);
	}
	getStatus();
};

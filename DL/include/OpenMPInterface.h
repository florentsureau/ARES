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
 **    File:  OpenMPInterface.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Interface to set OpenMP threads
 **    -----------
 ** - 09/2017: Added automatic/manual setting
 ******************************************************************************/

#ifndef OpenMPInterface_H
#define OpenMPInterface_H
#ifdef _OPENMP
	#include <omp.h>
#endif
#include <iostream>
#include "stdlib.h"
#include "math.h"


class OMPThreadHandling {

protected:
	static unsigned int mOuterThreads;//!< Number of outer loop threads
	static unsigned int mInnerThreads;//!< Number of inner loop threads
	static unsigned int mTotalThreads;//!< Total number of threads

	/**
      * Set the number of OpenMP threads to use.
       * @param[in] outer number of threads for outer loop.
	   * @param[in] inner number of threads for inner loop.
     */
	static void setNThreads(const int outer,const int inner);

public:
	//Constructor/Destructors
	~OMPThreadHandling(){};

	/**
      * Display OpenMP status
     */
	static void getStatus();
	/**
	* Automatically set the number of OpenMP threads to use.
	 * @param[in] outerMax: max desired number of outer threads.
	 * @param[in] innerMax: max desired number of inner threads.
	 * @param[in] nbProcs: desired number of processors to use.
     */
 	static void autoSetThreads(const int outerMax, const int innerMax,const int nbProcs=0);
	/**
      * Get number of inner threads
	   *@return number of inner threads
     */
	inline static unsigned int getNInnerThreads(){return mInnerThreads;}
	/**
	 * Get number of outer threads
	  *@return number of inner threads
     */
	inline static unsigned int getNOuterThreads(){return mOuterThreads;}

	/**
	  * Get total number of threads
	   *@return total number of threads
	 */
	inline static unsigned int getNThreads(){return mTotalThreads;}

};
#endif

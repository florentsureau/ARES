
#ifndef _IM_GLOB_H_
#define _IM_GLOB_H_

#include<cmath>
#include<cstdio>
// #include<cstring>
#include<cassert>
#include<cstdlib>
#include<iostream>
#include<string>
#include<sstream>
#include<fstream>

#include<complex>
#ifdef __WIN32__
#include <unistd.h>
#endif

//using namespace std;
typedef std::complex<float> complex_f;
typedef std::complex<double> complex_d;

// #include<climits>

#ifndef __WIN32__
#ifndef OSF1
#ifndef HP
#ifndef MACOS
#include <values.h> // deprecated, replaced by float & limits
//#include <float.h>
//#include <limits.h>
#endif
#endif
#endif
#endif

// On Windows machines, we need to define timersub.
#ifdef __WIN32__
#include <winsock.h>  //timeval on windows
inline void timersub(const timeval* tvp, const timeval* uvp, timeval* vvp)
{
  vvp->tv_sec = tvp->tv_sec - uvp->tv_sec;
  vvp->tv_usec = tvp->tv_usec - uvp->tv_usec;
  if (vvp->tv_usec < 0)
  {
	 --vvp->tv_sec;
	 vvp->tv_usec += 1000000;
  }
}
#endif


extern "C"
{
#include "fitsio2.h"
#undef True
#undef False
}

#define DEFBOOL 1
#ifdef DEFBOOL
#undef False
#undef True
enum Bool {False = 0,True = 1};
#else
#undef False
#undef True
#define True 1
#define False 0
#endif

// output for help  
#define OUTMAN stdout
#define WRITE_PARAM 0
#define MAX_NL 35000
#define MAX_NC 35000

inline void manline()
{
   fprintf(OUTMAN, "\n");
}

#if VMS
inline char *strdup(char *s1)
{
   int T = strlen(s1);
   char *Ret = new char[T];
   strcpy(Ret, s1);
   return(Ret);
}
#endif

#include "SoftInfo.h"
#include "OptMedian.h"
#include "DefMath.h"
#include "Memory.h"
#include "Array.h"
#include "Licence.h"
#include "Usage.h"

int GetOpt(int argc, char **argv, const char *opts);

#endif

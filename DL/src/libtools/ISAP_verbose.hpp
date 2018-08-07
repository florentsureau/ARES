/*============================================================================*/
/* File Description                                                           */
/*============================================================================*/
/**
 * @file        ISAP_verbose.hpp
 * @version     $Revision: 1.0 $
 * @author      A. WOISELLE
 *
 * Warning : these fuctions do not accept comas "," in the included code.
 * example :     " int a, b; "
 *  must become  " int a; int b; "
 *
 */
/*============================================================================*/

#ifndef verbose_macros_H_
#define verbose_macros_H_

#ifdef ENABLE_VERBOSE_ERROR
#define ERROR_(x) x
#else
#define ERROR_(x)
#endif

#ifdef ENABLE_VERBOSE_INFO
#define INFO_(x) x
#else
#define INFO_(x)
#endif

#ifdef ENABLE_VERBOSE_DEBUG
#define DEBUG_(x) x
#else
#define DEBUG_(x)
#endif

#define FINFO_ DEBUG_(std::cout << "> " << __PRETTY_FUNCTION__ << std::endl;)
#define FOFNI_ DEBUG_(std::cout << "< END " << __PRETTY_FUNCTION__ << std::endl;)

#endif // verbose_macros_H_

/*
 *
 * Authors: Christos Bakogiannis (mpakoyan "at" gmail "dot" com)
 *          Nikolaos Karapanos   (karapas "at" gmail "dot" com)
 *
 *
 * This file is part of kmGNFS.
 *
 * This distribution of the General Number Field Sieve algorithm is
 * placed in the public domain by its authors. It is distributed with
 * the hope that it will be useful, and you can use it freely for any
 * purpose without having to notify anyone. However, if you find our
 * work useful in any way, please let us know about it.
 *
 * kmGNFS makes use of Victor Shoup's NTL free library.
 * NTL can be used in conjunction with GMP (the GNU Multi-Precision library)
 * for enhanced performance. For more info, visit NTL's home page.
 *
 * We disclaim any responsibility for any bug.
 *
 */


#ifndef _GNFS_GLOBAL_H
#define	_GNFS_GLOBAL_H

#include <GNFS/gnfs_utils.h>


//paremeterized integer type that will be used by the algorithm in various places
typedef int int_type;
typedef unsigned int uint_type;
#define GNFS_BITS_PER_INT_TYPE NTL_BITS_PER_INT//NTL_BITS_PER_LONG //CAUTION: Change this as well!!!!!!



//here define some trully 64-bit types of the target platform
/*typedef long int64_type;
typedef unsigned long uint64_type;*///todo future not needed for the moment




//Various compile time parameters
#include <cassert>
#if 1
#define GNFS_DEBUG

/*
 * If GNGS_DEBUG is defined, then various informative messages will be printed and checks will be performed during the execution.
 * Bear in mind that, if some of these checks fail they will stop the execution of the program.
 * Disable this for improved performance.
 */

#endif


#if 1

#define GNFS_FACTOR_BASES_USE_FINDROOTS

/*
 * If FACTOR_BASES_USE_FINDROOTS is defined, then during the construction of AFB and QCB the FindRoots approach will be used.
 * If it is not defined the CanZass factoring algorithm will be used instead.
 */

#endif


#if 0

#define GNFS_USE_ORIGINAL_F1_F2

//ATTENTION: This flag is not used anymore! Only F2 = x - my and F1 = (y)^d*f(x/y) are used by this implementation!
/*
 * If GNFS_USE_ORIGINAL_F1_F2 is defined, then the polynomials F2 = x + my and F1 = (-y)^d*f(-x/y)  are used.
 * If it is not defined, the polynomials F2 = x - my and F1 = (y)^d*f(x/y) are used instead.
 */
//ATTENTION: This flag is not used anymore! Only F2 = x - my and F1 = (y)^d*f(x/y) are used by this implementation!

#endif


#if 1

#define GNFS_LINE_SIEVING_WITH_LOGARITHMS

/*
 * If GNFS_LINE_SIEVING_WITH_LOGARITHMS is defined, then line sieving with logarithms is used.
 * If it is not defined, simple line sieving with pure trial division is used instead.
 * WARNING!: line sieving with trial division is extremely SLOW and should be used only for testing purposes and for really small numbers!
 */

#endif

#endif	/* _GNFS_GLOBAL_H */


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


#ifndef _GNFS_UTILS_H
#define	_GNFS_UTILS_H

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/RR.h>

NTL_CLIENT

//Horner evaluation for ZZX class

void eval(ZZ& res, const ZZX& f, const ZZ& a); //CAUTION: res must not be an alias of a!!

inline ZZ eval(const ZZX& f, const ZZ& a) {

    ZZ x;
    eval(x, f, a);
    NTL_OPT_RETURN(ZZ, x);
}

inline void eval(ZZ& res, const ZZX& f, const long a) {

    static ZZ x;
    conv(x, a);
    eval(res, f, x);
}

inline ZZ eval(const ZZX& f, const long a) {

    static ZZ x;
    conv(x, a);
    return eval(f, x);
}


//Evaluation of the homogeneous polynomial F(x,y) = y^d * f(x/y) = c_d*x^d + c_d_1*x^(d-1)*y + ... + c_0*y

void evalHomogeneousF(ZZ& res, const ZZX& f, const ZZ& x, const ZZ& y); //CAUTION: res must not be an alias of x or y!!

inline void evalHomogeneousF(ZZ& res, const ZZX& f, const long x, const long y) {

    static ZZ a, b;
    conv(a, x);
    conv(b, y);
    evalHomogeneousF(res, f, a, b);

}

//Floating point version using doubles for input/output and RRs for the intermediate calculations
double evalHomogeneousF(const ZZX& f, const double x, const double y);

//Floating point version using only RRs
//CAUTION: res must not be an alias of x or y!!
void evalHomogeneousF(RR& res, const ZZX& f, const RR& x, const RR& y);


//print factorization of n
void printFactorizationOfN(const ZZ& n, const ZZ& aFactorOfN);


void printGNFSUsage();

//2D array dynamic allocation

template<class Item>
Item **alloc2D(long r, long c) {
    Item **t = new (nothrow) Item*[r];
    if (t == NULL)
        return NULL;

    //be careful: in case of error, previously allocated memory is not deallocated
    for (long i = 0; i < r; i++) {
        t[i] = new (nothrow) Item[c];
        if (t[i] == NULL)
            return NULL;
    }
    return t;
}


#endif	/* _GNFS_UTILS_H */


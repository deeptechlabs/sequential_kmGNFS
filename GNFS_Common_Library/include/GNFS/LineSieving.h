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


#ifndef _GNFS_LINE_SIEVING_H
#define	_GNFS_LINE_SIEVING_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <strings.h>
#include <sstream>
#include <string>
#include <signal.h>


#include <NTL/ZZ.h>
#include <NTL/ZZVec.h>

#include <GNFS/Sieving.h>
#include <GNFS/Relation.h>
#include <GNFS/gnfs_utils.h>
#include <GNFS/gnfs_global.h>

NTL_CLIENT


class LineSieveInfo : public SieveInfo {
public:

    int_type* b_x_r_mod_p;
    int_type* a_min_mod_p;

    sieve_type* row;
    long rowLength;

    LineSieveInfo() : SieveInfo(), b_x_r_mod_p(0), a_min_mod_p(0), row(0) {
    }

    ~LineSieveInfo() {
        if (b_x_r_mod_p)
            delete [] b_x_r_mod_p;
        if (a_min_mod_p)
            delete [] a_min_mod_p;
        if (row)
            delete [] row;
    }
};

class LineSieving : public Sieving {
public:

    static const char* relationsFileInfix;
    bool lineSevingInitialized;

    LineSieveInfo sieve_R;
    LineSieveInfo sieve_A;

    long last_b; //last b value that has been previously sieved

public:

    LineSieving()
    : Sieving(), lineSevingInitialized(false) {

    }

    bool init();

    bool doLineSieving() {
#ifdef GNFS_LINE_SIEVING_WITH_LOGARITHMS
        cout << "Applying line sieving with logarithms..." << endl;
        return doLineSievingWithLogarithms();
#else
        cout << "Applying line sieving with trial division..." << endl;
        return doSimpleLineSieving();
#endif
    }


protected:

    bool doLineSievingWithLogarithms();
    bool doSimpleLineSieving();
    bool trialDivideRelation(Relation& rel, const ZZ& b_x_m, const int_type* b_x_r_mod_p_R, const int_type* b_x_r_mod_p_A);
    bool writeRelationToFile(const Relation& rel);

    bool getSievingStatus();
    bool updateSievingStatus();


    long calculateSievingOffset(const int_type a_min_mod_p, const int_type b_x_m_mod_p, const long p);
    int_type reduce(const long a, const int_type p); //reduce a mod p

};

inline long LineSieving::calculateSievingOffset(const int_type a_mod_p, const int_type b_x_k_mod_p, const long p) {


    long sievingOffset = -SubMod(a_mod_p, b_x_k_mod_p, p);

#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING) && !defined(NTL_CLEAN_INT))
    sievingOffset += (sievingOffset >> (NTL_BITS_PER_LONG - 1)) & p;
    return sievingOffset;
#elif (defined(NTL_AVOID_BRANCHING))
    sievingOffset += (long) ((-(((unsigned long) sievingOffset) >> (NTL_BITS_PER_LONG - 1))) & ((unsigned long) p));
    return sievingOffset;
#else
    if (sievingOffset < 0)
        return sievingOffset + p;
    return sievingOffset;
#endif
}

inline int_type LineSieving::reduce(const long a, const int_type p) { // reduce a mod p

    int_type a_mod_p = a % p;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING) && !defined(NTL_CLEAN_INT))
    a_mod_p += (a_mod_p >> (GNFS_BITS_PER_INT_TYPE - 1)) & p;
    return a_mod_p;
#elif (defined(NTL_AVOID_BRANCHING))
    a_mod_p += (long) ((-(((unsigned long) a_mod_p) >> (GNFS_BITS_PER_INT_TYPE - 1))) & ((unsigned long) p));
    return a_mod_p;
#else
    if (a_mod_p < 0)
        return a_mod_p + p;
    return a_mod_p;
#endif

}

#endif	/* _GNFS_LINE_SIEVING_H */


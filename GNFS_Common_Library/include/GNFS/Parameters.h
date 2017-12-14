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


#ifndef _GNFS_PARAMETERS_H
#define	_GNFS_PARAMETERS_H

#include <iostream>
#include <fstream>
#include <climits>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include <GNFS/gnfs_global.h>

NTL_CLIENT

class Parameters {
public:
    ZZ n; //Number to factor
    long d; //Degree of polynomial
    bool monic; //true if the polynomial used by GNFS will be monic
    long maxRFB; //Rational factor base bound
    long maxAFB; //Algebraic factor base bound
    long QCBSize; //Quadratic character base size
    bool outputBasesInText; //Output bases in text format as well, if this is true.

    long max_a; //x dimension size (actually the half of it) of the line sieving step ( |a| <= value )
    long min_b; //Lower bound for y dimension of the line sieving step
    long max_b; //Upper bound for y dimension of the line sieving step (min_b <= b <= max_b)
    //int64_type L1; //Large prime bound for RFB (if 0, then the large prime optimization will not be used)//todo future
    //int64_type L2; //Large prime bound for AFB (if 0, then the large prime optimization will not be used)

    //ATTENTION!!: Large prime optimization not implemented yet!!
    long L1; //Large prime bound for RFB (if 0, then the large prime optimization will not be used) //todo future
    long L2; //Large prime bound for AFB (if 0, then the large prime optimization will not be used)
    long numOfExcessiveRelations; //Minimum number of excess relations to find during the sieving step

    ZZX f; //Polynomial f
    ZZX diffF; //Derivative of f

    ZZ m; //Root of f mod n
    ZZ minM; //Lower bound for m
    ZZ maxM; //Upper boumd for m


    /* Starting length in bits of the inert primes that will be used in the square root step.
     * The bigger this value is the fewer inert primes will be necessary (i.e. IPBsize can be smaller)
     * thus resulting in faster execution of the square root step.
     * Bear in mind that for every inert prime pi the following must hold: 
     * pi < NTL_SP_BOUND which also means that minInertPrimeLength <= NTL_SP_NBITS
     * (NTL_SP_BOUND is usually 2^30 on 32-bit machines and 2^50 on 64-bit machines)
     */
    long minInertPrimeLength;

    /* Number of inert primes that will be used for the square root step.
     * This value depends on  minInertPrimeLength and on n, i.e the number we are trying to factor.
     * If it is set too low then Couveignes square root algorithm will fail to work.
     * If it is set too high then the square root step will run longer that it is actually necessary.
     * Its value is best determined empirically.
     */
    long IPBSize;

    /* Maximum number of dependencies to process in the square root step.
     * A higher value means better odds to find a non trivial factorization of n
     * since there will be more dependencies to check but also more time to
     * compute the dependencies from the reduced matrix.
     * This should hold: maxNumOfDepsToProcess <= (num_of_cols_of_reduced_matrix - rank_of_reduced_matrix)
     * Each dependency has a 50% probability to yield a non trivial factorization
     * of n, so a reasonable value would be between 4 and 8.
     * If maxNumOfDepsToProcess <= 0
     * then the default of 4 dependencies will be used.
     * If maxNumOfDepsToProcess > (num_of_cols_of_reduced_matrix - rank_of_reduced_matrix)
     * then maxNumOfDepsToProcess = (num_of_cols_of_reduced_matrix - rank_of_reduced_matrix)
     */
    long maxNumOfDepsToProcess;

    static const char* paramsFile; //gnfs algorithm parameters

public:

    Parameters(bool readDefaultParamsFile = true) {

        if (readDefaultParamsFile)
            readInputParamsFile(paramsFile);
    }

    virtual ~Parameters() {
    }

    virtual void readInputFileStream(ifstream& istr);
    virtual void writeOutputFileStream(ofstream& ostr) const;
    void readInputParamsFile(const char* inputFile);
    void writeParamsFile() const;


};

#endif	/* _GNFS_PARAMETERS_H */


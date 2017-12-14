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


#ifndef _GNFS_INITIALIZATION_H
#define	_GNFS_INITIALIZATION_H

#include <iostream>
#include <fstream>
#include <climits>

#include <NTL/ZZ.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/RR.h>

#include <GNFS/Parameters.h>
#include <GNFS/gnfs_utils.h>

#include <GNFS/PolynomialSelection.h>

NTL_CLIENT


class Initialization : public Parameters {
public:
    bool extensivePolynomialSelection; //If true, extensive polynomial search will take place. If false, the algorithm will use the first polynomial that it will find, without any further search.
    bool calculateRFBAndAFBBound; //Programmatically calculate maxRFB and maxAFB
    bool performPolynomialSelection; //If true, perform the polynomial selection step. If false just skip it.


    double x1; //It is the lower bound (x1 <= |a_d/m|) defined in murhy's thesis for the non-monic polynomial selection
    double x2; //It is the upper bound (|a_d/m| <= x2) defined in murhy's thesis for the non-monic polynomial selection
    /*
     * Alternative way to handle the range of the leading coefficient's values.
     * If range > 0 then the x2 value is ignored and the search process will go on
     * from a_d_min up to a_d_min + range. It is used for faster polynomial searching,
     * but the resulted polynomial won't be the best
     */
    long range;

    /*
     * If false, the non-monic polynomial selection is performed by the least
     * alpha value. If true it is performed by the highest rating value.
     */
    bool byRating;

    static const char* inputParamsFile; //The initial parameters file

    long precision; //for RR class
    static const long maxAFBFactor = 2; // maxAFB = maxAFBFactor * maxRFB (if maxRFB and maxAFBL are calculated programmatically)
    static const double LExponent = 1.2; // L1 = maxRFB^LExponent, L2 = maxAFB^LExponent (if L1 and L2 are calculated programmatically)
    static const long numberOfMillerRabinTests = 100; //Miller-Rabin tests for the primality check of n.

    static const long verbose_ZZX_Factoring = 0; //1: NTL::factor executes in verbose mode, 0: NTL::factor does not execute in verbose mode
    static const long InitNumPrimes_ZZX_Factoring = 20; //Number of small primes p that are used in NTL::factor to factor f mod p (The NTL default is 7).


public:

    Initialization()
    : Parameters(false) {

        readInputParamsFile(inputParamsFile);
        initialCalculations();
    }

    void initialCalculations();

    void calculateMinM();

    void calculateMaxM();

    void calculateBounds();

    void calculateD();

    bool runInitialTests() const;

    bool isNPrime() const;

    bool testD() const;

    bool testM() const;

    bool polynomialSelectionStep();

    void readInputFileStream(ifstream& f);

    void writeOutputFileStream(ofstream& f) const;

};

#endif	/* _GNFS_INITIALIZATION_H */


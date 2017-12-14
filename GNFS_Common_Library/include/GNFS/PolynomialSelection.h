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


#ifndef _GNFS_POLYNOMIAL_H
#define	_GNFS_POLYNOMIAL_H

#include <vector>//todo
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/lzz_pXFactoring.h>

#include <queue>
#include <algorithm>

#include <GNFS/Parameters.h>

NTL_CLIENT

class ratedPolynomials {
    bool byRating;

public:

    ratedPolynomials(const bool& byRating_ = false) {
        byRating = byRating_;
    }


    ZZX poly;
    ZZ m;
    double alpha;
    double rating;

    bool operator() (const ratedPolynomials& lhs, const ratedPolynomials& rhs) const {
        if (byRating)
            return (lhs.rating > rhs.rating);
        else
            return (lhs.alpha < rhs.alpha);
    }
};

class ratedPolynomials_RR {
public:

    ratedPolynomials_RR() {
    }


    ZZX poly;
    ZZ m;
    RR alpha;
    RR rating;

    bool operator() (const ratedPolynomials_RR& lhs, const ratedPolynomials_RR& rhs) const {
        return (lhs.alpha < rhs.alpha);
    }
};

typedef priority_queue< ratedPolynomials, vector< ratedPolynomials >, ratedPolynomials > priorityQueueForRatingPolynomials;
typedef priority_queue< ratedPolynomials_RR, vector< ratedPolynomials_RR >, ratedPolynomials_RR > priorityQueueForRatingPolynomials_RR;

class PolynomialSelection {
public:

    ZZ m;
    ZZX selectedPolynomial;
    vector<ZZ> mValueofEachPolynomial;

    static const char* bestPolynomialsFile; //The file that stores the best polynomials found from the process.

private:
    ZZ n;
    long d;
    long precision; //for RR class //todo

public:

    PolynomialSelection(ZZ& n_, long d_) : n(n_), d(d_) {
    }

    ZZX findMonicPolynomial(const ZZ& m) const;

    void getPrimesLessThan100(vector<long>& primesLessThan100);

    //Methods with double precision
    double calcAlpha(ZZX &f, vector<long> &primes); //double precision
    double computeContP(ZZX& f, long p, ZZ& disc, int bound); //double precision
    double ratePolynomial(ZZX &f, double alpha, long B); //double precision
    double dickmanRho(double &u); //double precision

    //Methods with better precision
    RR calcAlpha_RR(ZZX &f, vector<long> &primes); //RR precision
    RR computeContP_RR(ZZX &f, RR &prime, ZZ &disc, int bound); //RR precision
    RR ratePolynomial_RR(ZZX &f, RR &alpha, long B); //RR precision
    RR dickmanRho_RR(RR &u); //RR precision

    struct checkRatingStruct {
        //When used in sort operation, this operator stores the polynomial with
        //the highest ranking in the first position of the vector (vector[0]).

        bool operator() (const ratedPolynomials& lhs, const ratedPolynomials & rhs) {
            return (lhs.rating > rhs.rating);
        }
    } checkRating;

    struct checkRatingStruct_RR {
        //When used in sort operation, this operator stores the polynomial with
        //the highest ranking in the first position of the vector (vector[0]).

        bool operator() (const ratedPolynomials_RR& lhs, const ratedPolynomials_RR & rhs) {
            return (lhs.rating > rhs.rating);
        }
    } checkRating_RR;


    bool findBestPolynomial(ZZX& f, ZZ& m, const ZZ& minM, const ZZ& maxM, const long B, const bool byRating, const long numOfBestRatedPolynomialsToStore, const bool performSieveTest);
    void findTopRatedPolynomials(vector<ratedPolynomials> &topRatedPolynomialsVector, const ZZ& minM, const ZZ& maxM, const long B, const bool byRating, long numOfBestRatedPolynomialsToStore);
    //void  sieveTestOnTopRatedPolynomials(pair<ZZX,ZZ> &topRatedPolynomials);

    bool findBestPolynomial_RR(ZZX& f, ZZ& m, const ZZ& minM, const ZZ& maxM, const long B, const long numOfBestRatedPolynomialsToStore, const bool performSieveTest);
    void findTopRatedPolynomials_RR(vector<ratedPolynomials_RR> &topRatedPolynomialsVector, const ZZ& minM, const ZZ& maxM, const long B, const long numOfBestRatedPolynomialsToStore);
    //void  sieveTestOnTopRatedPolynomials(pair<ZZX,ZZ> &topRatedPolynomials);

    ZZX findNonSkewedNonMonicPolynomial(const ZZ& m) const;
    void findBestNonSkewedNonMonicPolynomial(ZZX& f, ZZ& m, const long B, const ZZ& minM, const ZZ& maxM, const double x1, const double x2, const long range, const bool byRating, const long numOfBestRatedPolynomialsToStore, const bool performSieveTest);

    ZZX findNonSkewedNonMonicPolynomialWithLeadingCoefKnown(const ZZ& m, const ZZ &leadingCoef) const;
    void tweakRestCoeffsOfKnownNonSkewedNonMonicPolynomial(ZZX &poly, const ZZ& m, bool &changed);
    bool checkCoeffsSizeToDiscardPolynomial(ZZX &poly, const ZZ& m);


    //todo
    //ZZX tweakMonicPolynomial(Parameters& params, ZZX& coeffs, ZZ& m);
    //vector<ZZX> getAllPolynomialsZZX(Parameters& params);
    //ZZX testAllPolynomialsAndReturnTheBestOne(vector<ZZX>& allPolynomials);
    //void storeAllPolynomialsInFile(vector<ZZX>& allPolynomials, vector<ZZ>& mValueofEachPolynomial);

};



#endif	/* _GNFS_POLYNOMIAL_H */


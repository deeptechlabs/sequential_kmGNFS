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


#ifndef _GNFS_SIEVING_H
#define	_GNFS_SIEVING_H

#include <fstream>
#include <cmath>

#include <GNFS/gnfs_global.h>
#include <GNFS/FactorBases.h>

using namespace std;

//data type for sieving with logarithms. 
//uint8_t is probably the best choice
//CAUTION: Do not change cause the sieving program will break!! (check logTarget and most importantly memset operations in sieving)
typedef uint8_t sieve_type;

class SieveInfo {
public:

    vec_pair_intType_intType* FB;
    long FBSize;
    long whichBase;

    long startingPrimeOffset;

    long a_min;
    long a_max;
    long b_min;
    long b_max;


    double logBase;
    double logOfBase;
    sieve_type* logP;
    sieve_type projectiveRootsLogBias; //sum of logs of the projective roots that divide current b


    long fudgeFactor; //lamda * log(factor base bound);
    //long scaledFudgeFactor;

    SieveInfo() : logP(0), FB(0) {
    }

    virtual ~SieveInfo() {
        if (logP)
            delete [] logP;
    }
};

class Sieving : public FactorBases {
public:

    //fix this value such that the logarithm of the largest possible sieve value fits in a sieve_type
    //(in other words logTarget should be smaller that the maximum value that sieve_type can store
    static const long logTarget = 200;

    static const long logBaseUpdateRate = 100; //update log base every logBaseUpdateRate b (or y, whatever you like!)

    static const long cuttofUpdateRate_R = 10000;
    static const long cuttofUpdateRate_A = 1000; //cuttofUpdateRate_R % cuttofUpdateRate_A should be zero

    //do not sieve with primes below ignoreSmallPrimesBound
    //Sieving::init() will auto fix this bound if factor bases are too small
    //default value: 30
    static long ignoreSmallPrimesBound;

    static const double lambda = 1.5; //fudgeFactor = lamda * log(factor base bound); reasonable values: 1.0 - 2.6


    long numOfReports;
    long numOfRequiredRelations;
    long numOfSievingFoundRelations;
    long relationsFileCounter;

    static const long relationsPerFileBound = 100000;
    static const char* relationsFileDirectory;
    static const char* relationsFilePrefix;
    static const char* relationsFileSuffix;
    static const char* sievingStatusFile;

    string currentRelationsFilePath;
    string sievingStatusFilePath;
    string relationsFileDirectoryPath;
    ifstream relationFileInStream;
    ofstream relationFileOutStream;
    ifstream sievingStatusFileInStream;
    ofstream sievingStatusFileOutStream;

    bool sievingInitialized;

public:

    Sieving()
    : FactorBases(), sievingInitialized(false) {

    }

    //read bases from disk, check for previous sieving results etc
    virtual bool init(bool discardPreviousSievingResults = false);


protected:

    virtual bool getSievingStatus();
    virtual bool updateSievingStatus();
    virtual bool discardSievingStatus();

    long log(const long k, const double logOfBase) {
        //Custom base logarithm for longs. The result is rounded to the nearest integer.

        return static_cast<long> ((0.5 + (std::log(static_cast<double> (k)) / logOfBase)));
    }

    long log(const ZZ& k, const double logOfBase) {
        //This is the ZZ version for custom base logarithms. The result is rounded down.

        return static_cast<long> ((NTL::log(k) / logOfBase));
    }

    double getLogBase(const long a_min, const long a_max, const long b_max, const long whichBase);

    bool computeLogarithms(SieveInfo& sieve, const long b_max);

    //update one last time the sieving status file, and possibly other things as well.
    //call this from at the end of the main sieving functions
    virtual bool finalizeSieving();

};

#endif	/* _GNFS_SIEVING_H */


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


#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>

#include <NTL/ctools.h>

#include <GNFS/gnfs_global.h>
#include <GNFS/Sieving.h>




const char* Sieving::relationsFileDirectory = "IO/Sieving/Relations/";
const char* Sieving::relationsFilePrefix = "rels__";
const char* Sieving::relationsFileSuffix = ".txt";
const char* Sieving::sievingStatusFile = "_sieving_status_.txt";

long Sieving::ignoreSmallPrimesBound = 30;

NTL_CLIENT

bool Sieving::init(bool discardPreviousSievingResults) {

    assert(cuttofUpdateRate_R % cuttofUpdateRate_A == 0);


    cout << "Reading RFB from disk..." << endl;
    if (!readFBFromFile(FactorBases::RFB_)) {
        cout << "An error occured while reading RFB from disk..." << endl;
        return false;
    }

    cout << "Reading AFB from disk..." << endl;
    if (!readFBFromFile(FactorBases::AFB_)) {
        cout << "An error occured while reading AFB from disk..." << endl;
        return false;
    }

    cout << "Reading QCB from disk..." << endl;
    if (!readFBFromFile(FactorBases::QCB_)) {
        cout << "An error occured while reading QCB from disk..." << endl;
        return false;
    }


    if (!readNumberOfFreeRelations()) {
        cout << "An error occured while reading the free relations file..." << endl;
        return false;

    }


    if (NTL_OVERFLOW(1, RFBSize, AFBSize)) {
        cout << "RFBSize + AFBSize shoud be smaller than or equal to NTL_OVFBND! Aborting..." << endl;
        return false;
    }

    numOfRequiredRelations = RFBSize + AFBSize + 1; //one column for the sign
    if (!monic)
        ++numOfRequiredRelations; //one column which will consist only of 1's and will ensure that the number of dependent relations is even

    if (NTL_OVERFLOW(1, numOfRequiredRelations, QCBSize + numOfExcessiveRelations)) {
        cout << "RFBSize + AFBSize + QCBSize + numOfExcessRelations shoud be smaller than or equal to NTL_OVFBND! Aborting..." << endl;
        return false;
    }

    numOfRequiredRelations += QCBSize + numOfExcessiveRelations;

    cout << "Number of bases elements (rows of the matrix): " << numOfRequiredRelations - numOfExcessiveRelations << endl;
    cout << "Number of excessive relations to take into the matrix: " << numOfExcessiveRelations << endl;
    cout << "Total number of required relations (columns of the matrix): " << numOfRequiredRelations << endl;


    //fix paths of dirs and files - realativePath must have been set earlier if != "" 
    relationsFileDirectoryPath = relationsFileDirectory;
    sievingStatusFilePath = relationsFileDirectoryPath;
    sievingStatusFilePath += sievingStatusFile;

    if (discardPreviousSievingResults) {
        cout << "### Discarding previous sieving results... ###" << endl;
        if (!discardSievingStatus())
            return false;
    }
    else
        cout << "Checking for sieving results..." << endl;

    if (!getSievingStatus())
        return false;

    cout << "Number of free relations: " << numOfFreeRelations << endl;
    cout << "Number of found full relations (through sieving): " << numOfSievingFoundRelations << endl;
    cout << "Full relations file counter is " << relationsFileCounter << "..." << endl;
    cout << "Reports that have been previously tested: " << numOfReports << "..." << endl;


    //fix ignoreSmallPrimesBound if factor bases are small
    if (RFBSize <= 30)
        ignoreSmallPrimesBound = 0;
    else if (RFBSize <= 50)
        ignoreSmallPrimesBound = 6;
    else if (RFBSize <= 100)
        ignoreSmallPrimesBound = 12;
    else if (RFBSize <= 300)
        ignoreSmallPrimesBound = 20;

    sievingInitialized = true;

    return true;
}

bool Sieving::getSievingStatus() {

    numOfSievingFoundRelations = 0;
    relationsFileCounter = 0;
    numOfReports = 0;

    sievingStatusFileInStream.open(sievingStatusFilePath.c_str());
    if (!sievingStatusFileInStream.is_open()) {//we assume that the file does not exist, thus there are not any previous results
        sievingStatusFileInStream.clear();
        return true;
    }


    sievingStatusFileInStream.peek();
    if (sievingStatusFileInStream.eof()) {
        sievingStatusFileInStream.close();
        return true;
    }

    sievingStatusFileInStream >> numOfSievingFoundRelations;
    sievingStatusFileInStream.ignore(INT_MAX, '\n');
    sievingStatusFileInStream >> relationsFileCounter;
    sievingStatusFileInStream.ignore(INT_MAX, '\n');
    sievingStatusFileInStream >> numOfReports;
    sievingStatusFileInStream.ignore(INT_MAX, '\n');

    if (!sievingStatusFileInStream.good()) {
        cout << "Error while reading from " << sievingStatusFilePath << " ..." << endl;
        sievingStatusFileInStream.close();
        return false;
    }

    sievingStatusFileInStream.close();

    return true;
}

bool Sieving::updateSievingStatus() {

    if (!sievingStatusFileOutStream.is_open()) {

        sievingStatusFileOutStream.clear();
        sievingStatusFileOutStream.open(sievingStatusFilePath.c_str());
        if (!sievingStatusFileOutStream.is_open()) {
            cout << "Error opening " << sievingStatusFilePath << endl;
            cout << "Make sure that the directory " << relationsFileDirectoryPath << " exists and is accessible..." << endl;
            return false;
        }
    }
    else
        sievingStatusFileOutStream.seekp(0, ios::beg);

    sievingStatusFileOutStream << numOfSievingFoundRelations << "\t# Number of full relations found till now.\n";
    sievingStatusFileOutStream << relationsFileCounter << "\t# Number of full relation files.\n";
    sievingStatusFileOutStream << numOfReports << " #\tNumber of reports that were trial divided till now." << endl;

    if (!sievingStatusFileOutStream.good()) {
        cout << "Error while writing to " << sievingStatusFilePath << " ..." << endl;
        return false;
    }


    return true;
}

bool Sieving::discardSievingStatus() {

    DIR* pdir;
    struct dirent* entry;

    pdir = opendir(relationsFileDirectoryPath.c_str());

    if (!pdir) {
        cout << "Error opening directory " << relationsFileDirectoryPath << endl;
        cout << "Make sure that the directory exists and is accessible..." << endl;
        return false;
    }

    size_t prefixLength = strlen(relationsFilePrefix);

    errno = 0;
    while ((entry = readdir(pdir))) {
        if (strncmp(entry->d_name, relationsFilePrefix, prefixLength) == 0) {
            currentRelationsFilePath = relationsFileDirectoryPath;
            currentRelationsFilePath += entry->d_name;
            if (unlink(currentRelationsFilePath.c_str()) == -1) {
                cout << "Error while deleting \"" << currentRelationsFilePath << "\" ..." << endl;
                perror("The following error occured");
                closedir(pdir);
                return false;
            }

        }
    }
    closedir(pdir);

    if (errno) {
        cout << "Error while reading the contents of \"" << relationsFileDirectoryPath << "\" ..." << endl;
        perror("The following error occured");
        return false;
    }

    if (unlink(sievingStatusFilePath.c_str()) == -1) {
        if (errno != ENOENT) { //if the error is that the file does not exist, ignore it...
            cout << "Error while deleting \"" << sievingStatusFilePath << "\" ..." << endl;
            perror("The following error occured");
            return false;
        }
    }

    return true;
}

double Sieving::getLogBase(const long a_min, const long a_max, const long b_max, const long whichBase) {

    /* By Jason Papadopoulos:
       Decide on a base for the logs of one polynomial.

       The rational poly is assumed linear, its maximum value
       occurs at one of the endpoints of the sieve interval

       Rigorously finding the extreme values of the algebraic poly
       would require finding the minima and maxima, and comparing
       the polynomial values there to the values at a0 and a1. However,
       experiments show that for small b the values at the endpoints
       are much larger than those at the extreme values, and for large
       b the values are close but the extreme values tend to be outside
       the sieve interval. Hence we cheat and just do the same as with
       the rational poly */

    static ZZ b_x_m, F1, F2;
    double t;
    switch (whichBase) {

        case RFB_:
        {
            mul(b_x_m, m, b_max);


            sub(F1, a_min, b_x_m);
            abs(F1, F1);
            sub(F2, a_max, b_x_m);
            abs(F2, F2);

            if (F1 > F2)
                t = NumBits(F1);
            else
                t = NumBits(F2);

        }
            break;
        case AFB_:
        {


            evalHomogeneousF(F1, f, a_min, b_max);
            abs(F1, F1);
            evalHomogeneousF(F2, f, a_max, b_max);
            abs(F2, F2);
            /*static ZZ F3;
            evalHomogeneousF(F3, f, (a_max + a_min)/2, b_max);*/



            /*abs(F3, F3);
            if(F3 > F1 && F3 > F2)
                cout << F3 << ' '<<F1 << ' ' << F2<<endl;*/

            if (F1 > F2)
                t = NumBits(F1);
            else
                t = NumBits(F2);

        }
            break;
        default:
            return -1.0;

    }

    // the base to use is a number x such that log_x (2^t) = logTarget.
    return pow(2.0, t / logTarget);
}

bool Sieving::computeLogarithms(SieveInfo& sieve, const long b_max) {

    sieve.logBase = getLogBase(sieve.a_min, sieve.a_max, b_max, sieve.whichBase);
    sieve.logOfBase = std::log(sieve.logBase);

    vec_pair_intType_intType& FB = *(sieve.FB);

    long i;
    for (i = 0; i < sieve.FBSize; ++i) {
        sieve.logP[i] = static_cast<sieve_type> (log(FB[i].a, sieve.logOfBase));
    }

    sieve.fudgeFactor = static_cast<long> (lambda * log(FB[sieve.FBSize - 1].a, sieve.logOfBase));

    return true;
}

bool Sieving::finalizeSieving() {

    return updateSievingStatus();
}


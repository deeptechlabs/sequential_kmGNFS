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


#ifndef _GNFS_FACTORBASES_H
#define	_GNFS_FACTORBASES_H

#include <iostream>
#include <fstream>

#include <GNFS/Parameters.h>

#include <pair_intType_intType.h>
#include <vec_intType.h>

NTL_CLIENT

//The class name is a bit misleading because the class constructs all the required bases, like QCB for example, and not only the factor bases
class FactorBases : public Parameters {
public:

    enum Bases_ {
        RFB_ = 0, AFB_ = 1, QCB_ = 2
    };

    vec_pair_intType_intType RFB;
    vec_pair_intType_intType AFB;
    vec_pair_intType_intType QCB;

    long RFBSize;
    long AFBSize;

    long numOfFreeRelations;

    static const char* RFBFileText;
    static const char* RFBFileBinary;
    static const char* AFBFileText;
    static const char* AFBFileBinary;
    static const char* QCBFileText;
    static const char* QCBFileBinary;
    static const char* freeRelationsFile;

public:

    FactorBases()
    : Parameters() {

        numOfFreeRelations = -1;
    }

    bool createRFB();
    bool createAFB();
    bool createQCB();
    bool writeFBToFile(const long whichBase) const;
    bool readFBFromFile(const long whichBase);
    bool readNumberOfFreeRelations();


private:
    long binarySearch(const vec_pair_intType_intType& a, const int_type v, long l, long r);
    void sieveOfEratosthenes(vec_intType& sieve) const;
    bool writeFBBinary(ofstream& f, const vec_pair_intType_intType& FB) const;
    bool writeFBText(ofstream& f, const vec_pair_intType_intType& FB) const;
    bool readFBBinary(ifstream& f, vec_pair_intType_intType& FB);


};

#endif	/* _GNFS_FACTORBASES_H */


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


#include <NTL/vec_GF2.h>


#include <NTL/lzz_pEX.h>

#include <iostream>
#include <fstream>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/RR.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/pair_ZZX_long.h>

#include <GNFS/gnfs_utils.h>
#include <GNFS/gnfs_global.h>
#include <GNFS/FactorBases.h>
#include <GNFS/Relation.h>

#include <pair_intType_intType.h>
#include <vec_intType.h>
#include <triple.h>
#include <timer.h>


const char* FactorBases::RFBFileText = "IO/Bases/RFB_text.txt";
const char* FactorBases::RFBFileBinary = "IO/Bases/RFB_binary";
const char* FactorBases::AFBFileText = "IO/Bases/AFB_text.txt";
const char* FactorBases::AFBFileBinary = "IO/Bases/AFB_binary";
const char* FactorBases::QCBFileText = "IO/Bases/QCB_text.txt";
const char* FactorBases::QCBFileBinary = "IO/Bases/QCB_binary";
const char* FactorBases::freeRelationsFile = "IO/Bases/free_relations.txt";

NTL_CLIENT

long FactorBases::binarySearch(const vec_pair_intType_intType& a, const int_type v, long l, long r) {

    long m;
    while (r >= l) {
        m = (l + r) / 2;
        if (v == a[m].a)
            return m;
        if (v < a[m].a)
            r = m - 1;
        else
            l = m + 1;
    }
    return -1;
}

bool FactorBases::createRFB() {

    RR p_x = to_RR(maxRFB) / log(maxRFB); //approximation of the number of primes below maxRFB

    RFB.kill();
    RFB.SetMaxLength(to_long(1.1 * p_x));
    cout << "Estimation of RFB size: " << to_long(1.1 * p_x) << endl;

    vec_intType sieve;
    sieve.SetLength(maxRFB);

    sieveOfEratosthenes(sieve);

    pair_intType_intType pair;
    for (long i = 2; i < maxRFB; ++i) {//RFB consists of all primes up to maxRFB
        if (sieve[i]) {//i is prime
            pair.a = i;
            pair.b = m % i;
            append(RFB, pair);
        }
    }
    RFBSize = RFB.length();

    cout << "RFB consists of " << RFBSize << " pairs..." << endl;

    return true;
}

bool FactorBases::createAFB() {

    RR p_x = to_RR(maxAFB) / log(maxAFB); //approximation of the number of primes below maxAFB
    double avgRoots = d / 2.0; //crude average of the roots of f mod p. The real average is usually smaller
    AFB.kill();
    AFB.SetMaxLength(to_long(avgRoots * p_x));
    cout << "Estimation of AFB size: " << to_long(avgRoots * p_x) << endl;

    vec_intType sieve;
    sieve.SetLength(maxAFB);

    sieveOfEratosthenes(sieve);

    ofstream file(freeRelationsFile);
    if (!file.is_open()) {
        cout << "Error opening " << freeRelationsFile << "..." << endl;
        return false;
    }

    //timer t;
    //double finish;


    //used for storing free relations
    vec_Relation freeRelations;
    Relation rel;
    rel.relPair.a = 0;
    rel.relPair.b = 0;
    rel.sign = false;
    pair_intType_intType pair;
    triple_intType triple;

    numOfFreeRelations = 0;
    long length;


    //t.elapsed();

    for (long p = 2; p < maxAFB; ++p) { //AFB consists of all primes up to maxRFB
        if (sieve[p] == 0)//p is not a prime
            continue;

        if (divide(n, p)) {
            cout << "Impossible luck! We factored n while building AFB!" << endl;
            printFactorizationOfN(n, to_ZZ(p));
            exit(EXIT_SUCCESS);
        }

        zz_p::init(p, NextPowerOfTwo(d) + 1);

        zz_pX f_mod_p;
        conv(f_mod_p, f);

        //find and insert projective roots//
        if (coeff(f_mod_p, d) == 0) {
            pair.a = p;
            pair.b = p;
            append(AFB, pair);
        }

        MakeMonic(f_mod_p);

#ifdef GNFS_FACTOR_BASES_USE_FINDROOTS

        zz_pX X;
        SetX(X);

        zz_pX g;
        if (deg(f_mod_p) > 0)
            PowerXMod(g, p, f_mod_p);
        else
            SetCoeff(g, p, 1);

        sub(g, g, X);

        vec_zz_p roots;
        FindRoots(roots, GCD(f_mod_p, g));

        length = roots.length();
        for (long i = 0; i < length; ++i) {
            pair.a = p;
            pair.b = rep(roots[i]);
            append(AFB, pair);
        }

#else

        vec_pair_zz_pX_long factors;
        CanZass(factors, f_mod_p);
        //berlekamp(factors, f_mod_p);

        zz_p root;
        length = factors.length();
        for (long i = 0; i < length; ++i) {
            if (deg(factors[i].a) == 1) {
                NTL::negate(root, ConstTerm(factors[i].a));

                pair.a = p;
                pair.b = rep(root);
                append(AFB, pair);

            }
        }
#endif

        if (length == d && p < maxRFB && coeff(f_mod_p, d) != 0) {
            long AFBLength = AFB.length();

            rel.softReset();
            pair.a = p;
            pair.b = 1;

            long index = binarySearch(RFB, p, 0, RFBSize - 1);
            if (index == -1) {
                cout << "Could not find prime " << p << " in the rational factor base... This should never happen!" << endl;
                return false;
            }

            rel.appendRationalFactor(pair, index);

            triple.a = p;
            triple.c = 1;
            for (long i = AFBLength - d; i < AFBLength; ++i) {
                triple.b = AFB[i].b;
                rel.appendAlgebraicFactor(triple, i);
            }
            append(freeRelations, rel);
            ++numOfFreeRelations;
        }

        ///////////////////////////////////////////////////////////////////
    }

    /*finish = t.elapsed();
    cout << "Time for the CanZass method:  " << finish << endl;*/

    /*bool globalflag = true;
    if(AFB.length() != this->AFB.length())
    {
        cout << "OOOOPPSSSS!!" << endl;
        globalflag = false;
    }
    else
    {
        for(long i = 0 ; i < AFB.length() ; i++)
        {
            bool flag = false;
            for(long j = 0; j < AFB.length() ; j++)
            {
                if(this->AFB[i] == AFB[j])
                {
                    flag = true;
                    break;
                }
            }
            if(!flag)
            {
                cout << "OPPPPSSSS!!" << endl;
                globalflag = false;
                break;
            }
        }
        
        
    }
    if(globalflag)
        cout << "YYYYEEEESSSSSS!!!" << endl;*/


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (AFB[AFB.length() - 1].a >= NTL_SP_BOUND) {
        cout << "The prime numbers of the factor bases must be smaller than NTL_SP_BOUND (" << NTL_SP_BOUND << ")..." << endl;
        cout << "If the current system is 32-bit, please consider switching to a 64-bit one..." << endl;
        return false;
    }

    file << numOfFreeRelations << '\n';
    for (long i = 0; i < numOfFreeRelations; ++i) {
        file << freeRelations[i] << '\n';
    }

    if (!file.good()) {
        cout << "Error while writing to " << freeRelationsFile << "..." << endl;
        return false;
    }
    file.close();


    AFBSize = AFB.length();
    cout << "AFB consists of " << AFBSize << " pairs..." << endl;
    cout << "A total of " << numOfFreeRelations << " free relations were found..." << endl;

    return true;
}

bool FactorBases::createQCB() {

    QCB.kill();
    QCB.SetMaxLength(QCBSize);

    PrimeSeq seq;
    seq.reset(L2 + 1);

    //int64_type q;//todo future
    long q;
    long length;

    while (QCB.length() < QCBSize) {
        q = seq.next();
        if (q == 0) {
            cout << "Error while using PrimeSeq to generate the QCB. PrimeSeq cannot generate such large primes. Another method has to be implemented! Aborting.." << endl;
            exit(EXIT_FAILURE);
        }

        if (divide(coeff(f, d), q))
            continue;
        if (divide(n, q)) {
            cout << "Impossible luck! We factored n while bulding QCB!" << endl;
            printFactorizationOfN(n, to_ZZ(q));
            exit(EXIT_SUCCESS);
        }

        zz_p::init(q, NextPowerOfTwo(d) + 1);

        zz_pX f_mod_q;
        conv(f_mod_q, f);
        MakeMonic(f_mod_q);

        pair_intType_intType pair;


#ifdef GNFS_FACTOR_BASES_USE_FINDROOTS
        zz_pX X;
        SetX(X);
        zz_pX g = PowerXMod(q, f_mod_q) - X;

        vec_zz_p roots;
        FindRoots(roots, GCD(f_mod_q, g));

        length = roots.length();
        for (long i = 0; i < length; ++i) {
            if (QCB.length() >= QCBSize)
                break;

            if (divide(eval(diffF, rep(roots[i])), q))
                continue;

            pair.a = q;
            pair.b = rep(roots[i]);
            append(QCB, pair);

        }

#else

        vec_pair_zz_pX_long factors;
        CanZass(factors, f_mod_q);
        zz_p s;

        length = factors.length();
        for (long i = 0; i < length; ++i) {
            if (QCB.length() >= QCBSize)
                break;

            if (deg(factors[i].a) == 1) {
                NTL::negate(s, ConstTerm(factors[i].a));

                if (divide(eval(diffF, rep(s)), q))
                    continue;

                pair.a = q;
                pair.b = rep(s);
                append(QCB, pair);

            }
        }
#endif
    }

    if (QCB[QCB.length() - 1].a >= NTL_SP_BOUND) {
        cout << "The prime numbers of the QCB must be smaller than NTL_SP_BOUND (" << NTL_SP_BOUND << ")..." << endl;
        cout << "If the current system is 32-bit, consider switching to a 64-bit one..." << endl;
        return false;
    }

    QCBSize = QCB.length();

    cout << "QCB consists of " << QCBSize << " pairs..." << endl;



    return true;
}

bool FactorBases::writeFBToFile(const long whichBase) const {

    const char* outputBinaryFile;
    const char* outputTextFile;
    const vec_pair_intType_intType* FB;

    switch (whichBase) {
        case RFB_:
            FB = &RFB;
            outputBinaryFile = RFBFileBinary;
            outputTextFile = RFBFileText;
            break;
        case AFB_:
            FB = &AFB;
            outputBinaryFile = AFBFileBinary;
            outputTextFile = AFBFileText;
            break;
        case QCB_:
            FB = &QCB;
            outputBinaryFile = QCBFileBinary;
            outputTextFile = QCBFileText;
            break;
        default:
            cout << "Wrong base identifier..." << endl;
            return false;
    }

    ofstream f;
    f.open(outputBinaryFile, ios::binary);
    if (!f.is_open()) {
        cout << "Error opening " << outputBinaryFile << "..." << endl;
        return false;
    }

    cout << "Writing base in binary format..." << endl;
    if (!writeFBBinary(f, *FB)) {
        f.close();
        return false;
    }

    f.close();

    if (!outputBasesInText)
        return true;

    f.open(outputTextFile);

    if (!f.is_open()) {
        cout << "Error opening " << outputTextFile << "..." << endl;
        return false;
    }

    cout << "Writing base in text format as well..." << endl;
    if (!writeFBText(f, *FB)) {
        f.close();
        return false;
    }

    f.close();
    return true;
}

bool FactorBases::readFBFromFile(const long whichBase) {

    const char* inputBinaryFile;
    vec_pair_intType_intType* FB;

    switch (whichBase) {
        case RFB_:
            FB = &RFB;
            inputBinaryFile = RFBFileBinary;
            break;
        case AFB_:
            FB = &AFB;
            inputBinaryFile = AFBFileBinary;
            break;
        case QCB_:
            FB = &QCB;
            inputBinaryFile = QCBFileBinary;
            break;
        default:
            cout << "Wrong base identifier..." << endl;
            return false;
    }

    ifstream f;
    f.open(inputBinaryFile, ios::binary);
    if (!f.is_open()) {
        cout << "Error opening " << inputBinaryFile << "..." << endl;
        return false;
    }

    if (!readFBBinary(f, *FB)) {
        f.close();
        return false;
    }

    switch (whichBase) {
        case RFB_:
            RFBSize = FB->length();
            break;
        case AFB_:
            AFBSize = FB->length();
            break;
        case QCB_:
            QCBSize = FB->length();
            break;
        default:
            cout << "Wrong base identifier..." << endl;
            return false;
    }

    f.close();
    return true;

}

bool FactorBases::readNumberOfFreeRelations() {

    ifstream f;
    f.open(freeRelationsFile);
    if (!f.is_open()) {
        cout << "Error opening " << freeRelationsFile << "..." << endl;
        return false;
    }

    f >> numOfFreeRelations;

    if (!f.good()) {
        cout << "Error while reading from " << freeRelationsFile << "..." << endl;
        f.close();
        return false;
    }

    f.close();
    return true;
}

void FactorBases::sieveOfEratosthenes(vec_intType& sieve) const {

    long i;
    long length = sieve.length();
    for (i = 2; i < length; ++i)
        sieve[i] = 1;
    for (i = 2; i * i < length; ++i) {
        if (sieve[i]) {
            for (long j = i * i; j < length; j += i)
                sieve[j] = 0;
        }
    }

}

bool FactorBases::writeFBBinary(ofstream& f, const vec_pair_intType_intType& FB) const {

    f.clear();
    long length = FB.length();
    f.write(reinterpret_cast<char*> (&length), sizeof (long));
    for (long i = 0; i < length; ++i) {
        f.write(reinterpret_cast<const char*> (&(FB[i].a)), sizeof (int_type));
        f.write(reinterpret_cast<const char*> (&(FB[i].b)), sizeof (int_type));
    }

    if (f.good())
        return true;
    return false;

}

bool FactorBases::writeFBText(ofstream& f, const vec_pair_intType_intType& FB) const {

    f.clear();
    long length = FB.length();
    f << length << endl;
    for (long i = 0; i < length; ++i) {
        f << FB[i] << endl;
    }

    if (f.good())
        return true;
    return false;
}

bool FactorBases::readFBBinary(ifstream& f, vec_pair_intType_intType& FB) {

    f.clear();
    long length;
    f.read(reinterpret_cast<char*> (&length), sizeof (long));

    //vec_pair_intType_intType RFB;
    FB.kill();
    FB.SetLength(length);

    for (long i = 0; i < length; ++i) {
        f.read(reinterpret_cast<char*> (&(FB[i].a)), sizeof (int_type));
        f.read(reinterpret_cast<char*> (&(FB[i].b)), sizeof (int_type));
    }

    /*if(this->QCB == FB )
        cout << "OOOOKKKKKKK!" << endl;
    cout << this->QCB << endl;
    cout << FB << endl;*/

    if (f.good())
        return true;
    return false;

}

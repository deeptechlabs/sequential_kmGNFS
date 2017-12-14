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


#ifndef _GNFS_SQUAREROOT_H
#define	_GNFS_SQUAREROOT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <strings.h>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pXFactoring.h>

#include <GNFS/gnfs_global.h>
#include <GNFS/LinearAlgebra.h>
#include <GNFS/Relation.h>



#if defined(GNFS_DEBUG) && 0
#define GNFS_SQRT_STEP_DEBUG

/*
 * If GNFS_ROOT_STEP_DEBUG is defined, then verbose informative messages will be printed and checks will be performed
 * during the execution of the square root step.
 * Bear in mind that, if some of these checks fail they will stop the execution of the program.
 *
 * CAUTION: Enabling this flag may cause severe performance penalty!
 * Use this only for debugging the square root step code using fairly small instances of the problem!
 */

#endif

class SquareRoot : public LinearAlgebra {
public:
    vec_vec_Relation dependencies;
    mat_GF2 solutions;

    //g(x) = f(x/c_d)*(c_d^(d - 1)), where c_d is the leading coefficient of f
    //g(x) defines the number ring Z[w]/g(w) that we will work into (w = c_d * a, a is a complex root of f)
    ZZX g;
    ZZ discr_g;

    vector<long> IPB; //primes p1,p2,...,pk that will be used as finite fields (Inert Prime Base)
    vector<long> xi_vec; //corresponding xis of the primes that will be used in the chinese remainder theorem
    vector<long> ai_vec; //ai = Pi^(-1) mod pi, corresponding to the above primes pi (Pi = P/pi, P = p1*p2*..*pk)

    long prime; //currently processed prime
    long primeIndex; //index in vector "primes" of currently processed prime (it will be used to properly store and retrieve the corresponding xi)
    ZZ q; //prime^d
    ZZ q_minus_1; //q - 1

    //map<long, zz_pE> xis;
    unsigned short *factor_exps_R; //exponents of the factors on the rational side
    unsigned short *factor_exps_A; //exponents of the factors on the algebraic side
    ZZ squareRootOfNormOfDep_A; //the square root in Z of the algebraic norm of a dependency

    bool squareRootInitialized;

    static const char* outputFile;

public:

    SquareRoot()
    : LinearAlgebra(), squareRootInitialized(false), factor_exps_R(0), factor_exps_A(0) {

    }

    ~SquareRoot() {
        if (factor_exps_R)
            delete [] factor_exps_R;
        if (factor_exps_A)
            delete [] factor_exps_A;
    }

    bool init(bool reducedMatrixBinaryInput = true);

    bool executeSquareRootStep();


private:

    long determineApplicableFiniteFields();
    void calculateSquareRootOfNormOfDep_A(const vec_Relation& dependency);
    bool findSquareRootOfFiniteField(const vec_Relation& dependency, const zz_pX& f_mod_p, const zz_pX& g_mod_p);
    bool findPolynomial_mod_f_p(zz_pE& dp, const vec_Relation& dependency, const zz_pX& f_mod_p, const zz_pX& g_mod_p); //Calculates Dp [f'(a)*PROD(a+ba)(mod p_i)] or else Eq. 1 Couveignes modulo prime
    void calculate_r_s(long& r, ZZ& s);
    bool findQuadraticNonResidue(zz_pE& QNR);
    bool calculateSquareRootOfFiniteField(zz_pE& root, const zz_pE& dp);
    void calculateNormBp(zz_p& norm_mod_p, const zz_pE& squareRoot);
    void calculateNormBModPrime(zz_p& norm_mod_p, const zz_pX& f_mod_p, const long numOfElementsInDependency);
    void evaluate_xi(zz_p& xi, const zz_pE& squareRoot, const zz_pX& f_mod_p);

    void crt(ZZ_p& x_mod_n, const ZZ_p& P_mod_n);

    //void calculate_P(ZZ& P); useless
    void calculate_P_mod_n(ZZ_p& P_mod_n);

    void calculateVectorOf_Pi_mod_p(vector<long>& Pi_mod_p_vec);
    long calculate_Pi_mod_p(const long prime);

    void calculateVectorOf_ai(const vector<long>& Pi_mod_p_vec);
    long calculate_ai(const long Pi, const long prime);

    void calculate_r_mod_n(ZZ_p& r_mod_n);
    void calculate_x_mod_n(ZZ_p& x_mod_n, const ZZ_p& r_mod_n, const ZZ_p& P_mod_n);
    void calculateSumOfProduct_ai_xi_Pi_mod_n(ZZ_p& sumOfProduct_mod_n, const ZZ_p& P_mod_n);
    void calculateProductOf_r_P_mod_n(ZZ_p& r_x_P_mod_n, const ZZ_p& r_mod_n, const ZZ_p& P_mod_n);

    void findSquareRootInZ_mod_n(ZZ_p& squareRoot_mod_n, const vec_Relation& dependency);
    void calculateSquareRootOfDepInZ_mod_n(ZZ_p& squareRoot_mod_n, const vec_Relation& dependency);
    void calculate_diff_f_mod_n_in_m(ZZ_p& derivSqrt_mod_n, const ZZ_pX& f_mod_n);

    void calculate_gx(ZZX& g); //calculates g(x) = f(x/c_d)*(c_d^(d - 1)), where c_d is the leading coefficient of f

    long checkResult(ofstream& ofstr, const long currentDep, const ZZ& factor1, const ZZ& factor2);
    void outputFailureMsg1(ostream& ostr);
    void outputFailureMsg2(ostream& ostr);
};

#endif	/* _GNFS_SQUAREROOT_H */


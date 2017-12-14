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


#include <GNFS/SquareRoot.h>

#include <timer.h>

const char* SquareRoot::outputFile = "IO/SquareRoot/results.txt";

NTL_CLIENT

bool SquareRoot::init(bool reducedMatrixBinaryInput) {

    if (!LinearAlgebra::init())
        return false;

    cout << "Reading reduced matrix from disk..." << endl;
    if (reducedMatrixBinaryInput) {
        if (!readMatrixBinary(reducedMatrixBinaryFile))
            return false;
    }
    else {
        if (!readMatrixText(reducedMatrixTextFile))
            return false;
    }

    maxNumOfDepsToProcess = min(maxNumOfDepsToProcess, gnfsMatrix.NumCols() - gnfsMatrixRank);

    cout << "Matrix dimensions: " << gnfsMatrix.NumRows() << " x " << gnfsMatrix.NumCols() << endl;
    cout << "Matrix rank: " << gnfsMatrixRank << endl;
    cout << "Number of solutions to fetch: " << maxNumOfDepsToProcess << endl;

    cout << "Getting solutions and dependent relations..." << endl;

    if (!getSolutions(solutions, maxNumOfDepsToProcess))
        return false;

    if (!getDependenciesFromSolutions(dependencies, solutions))
        return false;

#ifdef GNFS_SQRT_STEP_DEBUG
    //this piece of code is a mesh up. it was used only for debugging in the early stages of development

    transpose(gnfsMatrix, gnfsMatrix);
    int_type *factors = new int_type[AFBSize];

    cout << "Checking fetched solutions for correnctness..." << endl;
    cout << "%%% ATTENTION!: This could take forever when dealing with big numbers! Use only for debuging small instances! %%%" << endl;
    for (long n = 0; n < solutions.NumRows(); ++n) {

        const vec_GF2& solution = solutions[n];
        //cout << endl << "solution:   " << solution << endl;
        //cout << solution << endl;
        vec_GF2 sum;
        sum.SetLength(gnfsMatrix.NumCols());

        for (long i = 0; i < solution.length(); ++i) {
            if (solution[i] == 1)
                sum += gnfsMatrix[i];
        }

        assert(weight(sum) == 0);



        const vec_Relation& dependency = dependencies[n];

        ZZ norm;
        ZZ tmp1;
        set(tmp1);
        ZZ tmp2;
        set(tmp2);
        vec_GF2 qcbcheck;
        ZZ tmpZZ, q;
        qcbcheck.SetLength(QCBSize);

        memset(factors, 0, AFBSize * sizeof (int_type));

        for (long i = 0; i < dependency.length(); ++i) {

            // cout << dependency[i].relPair << endl;
            const Relation& rel = dependency[i];
            long a = rel.relPair.a;
            long b = rel.relPair.b;


            if (b == 0)
                mul(tmp1, tmp1, rel.rationalFactors[0].a);
            else
                mul(tmp1, tmp1, a - b * m);


            if (b == 0)
                power(norm, rel.rationalFactors[0].a, d);
            else
                evalHomogeneousF(norm, f, a, b);
            mul(tmp2, tmp2, norm);
            //mul(tmp2, tmp2, evalHomogeneousF(f, a, b));
            //mul(tmp2, tmp2, evalHomogeneousF(f, a, -b));


            if (rel.relPair.b == 0) {//handle free relations separately
                long p = rel.rationalFactors[0].a; //get corresponing p

                for (long j = 0; j < QCBSize; ++j) {

                    conv(tmpZZ, p);
                    conv(q, QCB[j].a);
                    rem(tmpZZ, tmpZZ, q); //reduce numerator modulo q

                    if (Jacobi(tmpZZ, q) == -1L)// check the Legendre symbol of (p / q)
                        qcbcheck[j]++;

                }
            }
            else {
                for (long j = 0; j < QCBSize; ++j) {
                    conv(tmpZZ, QCB[j].b); // s
                    mul(tmpZZ, tmpZZ, b); // bs

                    sub(tmpZZ, a, tmpZZ); // a - bs

                    conv(q, QCB[j].a);
                    rem(tmpZZ, tmpZZ, q); //reduce numerator modulo q

                    /* if (Jacobi(tmpZZ, q) == 0) {
                         cout << "Jacobi: " << Jacobi(tmpZZ, q) << endl;
                         cout << tmpZZ << "  " << QCB[j].b << "  " << q << endl;

                     }*/

                    if (Jacobi(tmpZZ, q) == -1L)
                        qcbcheck[j]++; // = qcbcheck[j] + 1;

                }
            }

            for (long k = 0; k < rel.algebraicFactors.length(); ++k) {
                factors[rel.algebraicFactorsIndex[k]] += rel.algebraicFactors[k].c;
            }


            /*   ZZ NormAlg2;
               set(NormAlg2);
               ZZ tmp;
               for (long k = 0; k < AFBSize; ++k) {
                   power(tmp, AFB[k].a, factors[k]);
                   mul(NormAlg2, NormAlg2, tmp);
               }

               if (b == 0)
                  power(norm, rel.rationalFactors[0].a, d);
              else
                  evalHomogeneousF(norm, f, a, b);
            
               if(NormAlg2 != abs(norm)) {
                   cout << "ooopps. " << rel.relPair.a << ' ' << rel.relPair.b << endl;
                   cout << norm << ' ' << NormAlg2 << endl;
                    exit(1);
               }

               memset(factors, 0, AFBSize * sizeof (int_type));*/
        }
        /*cout << "[";
        for (long k = 0; k < AFBSize; ++k) {
            cout << factors[k] << ' ';
        }
        cout << "]" << endl;*/

        for (long k = 0; k < AFBSize; ++k) {
            if (factors[k] % 2 != 0)
                cout << "opps index: " << k << " ,prime: " << AFB[k] << endl;
        }

        ZZ rootNormAlg2;
        set(rootNormAlg2);
        ZZ tmp;
        for (long k = 0; k < AFBSize; ++k) {
            power(tmp, AFB[k].a, factors[k]);
            mul(rootNormAlg2, rootNormAlg2, tmp);
        }
        assert(rootNormAlg2 == tmp2);

        //cout <<rootNormAlg2 << endl;
        //cout << tmp2 << endl;
        //tmp2 = rootNormAlg2;

        ZZ rootZ = SqrRoot(tmp1);

        assert(rootZ * rootZ == tmp1);

        ZZ rootNormAlg = SqrRoot(tmp2);

        assert(rootNormAlg * rootNormAlg == tmp2);

        assert(weight(qcbcheck) == 0);

        // exit(1);
    }
    delete [] factors;
#endif

    gnfsMatrix.kill(); //the matrix is of no use anymore so give the memory back to the system.

    if (factor_exps_R == 0) {
        factor_exps_R = new (nothrow) unsigned short[RFBSize];
        if (factor_exps_R == 0) {
            cout << "Memory allocation error!!" << endl;
            return false;
        }
    }

    if (factor_exps_A == 0) {
        factor_exps_A = new (nothrow) unsigned short[AFBSize];
        if (factor_exps_A == 0) {
            cout << "Memory allocation error!!" << endl;
            return false;
        }
    }

    calculate_gx(g);
    discriminant(discr_g, g, 1);

    squareRootInitialized = true;
    return true;
}

bool SquareRoot::executeSquareRootStep() {


    ofstream ofstr;
    ofstr.open(outputFile);
    if (!ofstr.is_open()) {
        cout << "Could not open output file \"" << outputFile << "\" ..." << endl;
        return false;
    }
    long res;
    //cout << "f -> " << f << endl;

    /*
     * The first step is to find the primes that define applicable finite fields,
     * i.e. primes p_i for which f is irreducible modulo p_i.
     */
    cout << "Searching for inert primes pi to use..." << endl;
    res = determineApplicableFiniteFields();
    switch (res) {
        case -1 :
            return false;
        case 0:
            break;
        default:
            ZZ factor1, factor2;
            conv(factor1, res);
            div(factor2, n, factor1);
            checkResult(ofstr, -1, factor1, factor2);
            return true;
    }

    ZZ_p::init(n);


    //We calculate "P" as described in Couveignes algorithm
    /*cout << "Computing P..." << endl; //useless
    ZZ P;
    calculate_P(P);*/
    ofstr << "Square root step for n = " << n << " started ..." << endl << "------------------------------" << endl;

    cout << "Computing P mod n ( P = p1*p2*...*pk) ..." << endl;
    ZZ_p P_mod_n;
    calculate_P_mod_n(P_mod_n);
    //cout << "P -> " << P << endl;


    cout << "Computing Pi mod pi (Pi = P/pi) for all primes pi that will be used..." << endl;
    //We calculate "Pi"as described in Couveignes algorithm
    vector<long> Pi_mod_p_vec;
    calculateVectorOf_Pi_mod_p(Pi_mod_p_vec);



    cout << "Computing ai = Pi^(-1) mod pi for all primes pi that will be used..." << endl;
    //We calculate "ai"as described in Couveignes algorithm
    calculateVectorOf_ai(Pi_mod_p_vec);

    Pi_mod_p_vec.clear(); //not needed anymore


    long vecLength = IPBSize;
    xi_vec.resize(vecLength);

    long maxRoot = NextPowerOfTwo(d) + 1; //for zz_p::init


    for (long idx = 0; idx < dependencies.length(); ++idx) {
        ofstr << "Working with dependency " << idx + 1 << " of " << dependencies.length() << " ..." << endl;
        cout << "Working with dependency " << idx + 1 << " of " << dependencies.length() << " ..." << endl;


        const vec_Relation& dependency = dependencies[idx];

        calculateSquareRootOfNormOfDep_A(dependency);


        /*
         * Now, for each of the selected primes, we have to compute the square root in
         * the corresponding finite field (Fp^d), defined by the prime.
         */

        //timer ti;
        //ti.elapsed();

        for (long i = 0; i < vecLength;) {
            primeIndex = i;
            prime = IPB[i++]; //We next define a finite field using this prime
            power(q, prime, d);
            sub(q_minus_1, q, 1);

            if (i % 100 == 0) {
                cout << "Processing prime " << i << " of " << vecLength << endl;
            }


#ifdef GNFS_SQRT_STEP_DEBUG
            cout << endl << "################# iteration for prime " << prime << " #################" << endl;
#endif

            zz_p::init(prime, maxRoot);
            zz_pX f_mod_p;
            conv(f_mod_p, f);
            zz_pX g_mod_p;
            conv(g_mod_p, g);

            zz_pE::init(g_mod_p);


            if (!findSquareRootOfFiniteField(dependency, f_mod_p, g_mod_p))
                return false;


#ifdef GNFS_SQRT_STEP_DEBUG
            cout << "#############################################################" << endl;
#endif
        }
        //double t1 = ti.elapsed();
        //cout << "specific time: " << t1 << endl;

        /*
         * We can now use the Chinese Remainder Theorem to compute  b = B(a) e Z[a]
         * and \phi(\beta)= (B(m) mod n)e Z/nZ.
         */
        //timer ti;
        //ti.elapsed();

        ZZ_p x_mod_n;
        crt(x_mod_n, P_mod_n);

        //double t1 = ti.elapsed();
        //cout << "specific time: " << t1 << endl;

        cout << "x_mod_n -> " << x_mod_n << endl;
        ofstr << "x_mod_n -> " << x_mod_n << endl;



        //We now compute the square root in Z
        ZZ_p squareRootInZ_mod_n;
        findSquareRootInZ_mod_n(squareRootInZ_mod_n, dependency);


        cout << "squareRootInZ mod n -> " << squareRootInZ_mod_n << endl;
        ofstr << "squareRootInZ mod n -> " << squareRootInZ_mod_n << endl;

        cout << "Computing GCDs..." << endl;
        ZZ factor1, factor2;
        GCD(factor1, rep(x_mod_n) - rep(squareRootInZ_mod_n), n);
        GCD(factor2, rep(x_mod_n) + rep(squareRootInZ_mod_n), n);

        res = checkResult(ofstr, idx, factor1, factor2);
        switch (res) {
            case -1 :
                return false;
            case 0:
                continue;
            case 1:
                return true;
        }
    }

    return false;
}

long SquareRoot::determineApplicableFiniteFields() {

    //Populate the Inert Prime Base, that consists of primes p where polynomial f remains irreducible modulo p

    /* We could move tis code to the FactorBases class, so that the IPB is created only once together with the factor bases
     * and just read it from a file at this point, however the IPB creation is not that time consuming and poses no serious
     * problems the way it is implemented now. Maybe we 'll move the code in the future. //todo future
     */

    cout << "Need to find a total of " << IPBSize << " inert primes..." << endl;
    cout << "The starting length of the primes is " << minInertPrimeLength << " bits..." << endl;


    long maxRoot = NextPowerOfTwo(d) + 1; //for zz_p::init()

    IPB.clear();
    IPB.reserve(IPBSize);

    long p = 1L << (minInertPrimeLength - 1); //starting value for the search of inert primes
    p--;

    // timer ti;
    // ti.elapsed();
    for (long i = 0; i < IPBSize;) {

        p += 2;
        p = NextPrime(p, 50); //50 Miller-Rabin tests mean that the odds of a composite number passing the tests is 2^(-100)

        zz_p::init(p, maxRoot);
        zz_pX f_mod_p;
        conv(f_mod_p, f);


        //if (!IterIrredTest(f_mod_p))
        if (!DetIrredTest(f_mod_p))
            continue;

        if (f.rep[d] % p == 0)
            continue;


        if (n % p == 0 && n != p) {
            cout << "Great luck!! Prime " << p << " is a factor of n!" << endl;
            //printFactorizationOfN(n, to_ZZ(p));
            return p;
        }

        IPB.push_back(p);
        ++i;
    }

    // double tim = ti.elapsed();
    // cout << "specific time: " << tim << endl;

    assert(IPB.size() == IPBSize);
    return 0;
}

void SquareRoot::calculateSquareRootOfNormOfDep_A(const vec_Relation& dependency) {

#ifdef GNFS_SQRT_STEP_DEBUG
    ZZ prod, tmp, squareRoot;
    set(prod);
#endif

    bzero(factor_exps_A, AFBSize * sizeof (unsigned short));
    for (long i = 0; i < dependency.length(); ++i) {
        const Relation& rel = dependency[i];
#ifdef GNFS_SQRT_STEP_DEBUG
        long a = dependency[i].relPair.a;
        long b = dependency[i].relPair.b;
        evalHomogeneousF(tmp, f, a, b);
        prod *= tmp;
#endif

        long factorsLength = rel.algebraicFactors.length();
        for (long k = 0; k < factorsLength; ++k) {
            factor_exps_A[rel.algebraicFactorsIndex[k]] += rel.algebraicFactors[k].c;
        }
    }

    ZZ tmpCoFactor;
    set(squareRootOfNormOfDep_A);

    for (long k = 0; k < AFBSize; ++k) {
        power(tmpCoFactor, AFB[k].a, factor_exps_A[k] / 2);
        mul(squareRootOfNormOfDep_A, squareRootOfNormOfDep_A, tmpCoFactor);
    }


#ifdef GNFS_SQRT_STEP_DEBUG
    SqrRoot(squareRoot, prod);
    assert(squareRootOfNormOfDep_A == squareRoot);
#endif
}

bool SquareRoot::findSquareRootOfFiniteField(const vec_Relation& dependency, const zz_pX& f_mod_p, const zz_pX& g_mod_p) {


#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "=> In findSquareRootOfFiniteField()" << endl;
#endif

    //first compute the polynomial d_p (Briggs thesis)[or g (gamma - Couveignes)]
    //in the Finite Field defined by the prime.
    zz_pE dp;
    if (!findPolynomial_mod_f_p(dp, dependency, f_mod_p, g_mod_p))
        return false;

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "Dp -> " << dp << endl;
#endif



    //calculate the square root of dp in the finite field defined by the prime
    zz_pE root;
    if (!calculateSquareRootOfFiniteField(root, dp))
        return false;


#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "The first root is -> " << root << endl;
#endif

    /*
     * Next we have to make sure that the various roots bi that we calculate
     * (where bi = +-b mod prime) are congruent to the same square root b modulo p_i.
     * To do this we calculate the norm of bi (N(bi) = bi^((p_i^d - 1)/(p - 1) and the
     * residue modulo p_i of the norm of b  (N(b) as defined in Eq.(2) in Couveignes),
     * and compare the two norms. If the two norms are equal then bi=b mod p_i and we
     * keep Bi. If they are opposite, then we replace Bi by -Bi
     */

    //calculate the norm of the square root of dp we just found
    zz_p normOfRoot;
    calculateNormBp(normOfRoot, root);
#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "The norm of this root is -> " << normOfRoot << endl;
#endif

    //We calculate the residue modulo prime of the norm of b [N(a), Eq. (2) - Couveignes)
    zz_p normN_a;
    calculateNormBModPrime(normN_a, f_mod_p, dependency.length());
#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "N(a) -> " << normN_a << endl;
#endif


    zz_p xi;
    /*
     * iff the two norms are equal then the square root "root" is the wanted, otherwise
     * we need the "-root"
     */

    if (normN_a == normOfRoot) {
        evaluate_xi(xi, root, f_mod_p);
#ifdef GNFS_SQRT_STEP_DEBUG
        cout << "We calculate the xi for \"root\"" << endl;
        cout << "The xi for \"root\" is -> " << xi << endl;
#endif
    }
    else {
        NTL::negate(root, root); //We take "-root"



#ifdef GNFS_SQRT_STEP_DEBUG
        cout << "The -root is -> " << root << endl;
        calculateNormBp(normOfRoot, root);
        assert(normOfRoot == normN_a);
        cout << "The norm of  -root is -> " << normOfRoot << endl;
#endif

        evaluate_xi(xi, root, f_mod_p);
#ifdef GNFS_SQRT_STEP_DEBUG
        cout << "We calculate the xi for \"-root\"" << endl;
        cout << "The xi for \"-root\" is -> " << xi << endl;
#endif
    }

    xi_vec[primeIndex] = xi._zz_p__rep; //store the xi of the current prime in the corresponding cell

    return true;
}

bool SquareRoot::findPolynomial_mod_f_p(zz_pE& dp, const vec_Relation& dependency, const zz_pX& f_mod_p, const zz_pX& g_mod_p) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "==> In findPolynomial_mod_f_p()" << endl;
#endif

    // zz_pX dpp;
    // set(dpp);

    zz_p cd_x_a, a_mod_p, b_mod_p;
    // zz_pX smoothPair_mod_p;
    zz_pE smoothPair_mod_f_p;

    set(dp); //dp = 1

    long length = dependency.length();
    for (long i = 0; i < length; ++i) {
        //cout << dependency[i].relPair << endl;
        conv(a_mod_p, dependency[i].relPair.a);
        conv(b_mod_p, -dependency[i].relPair.b);

        mul(cd_x_a, f_mod_p.rep[d], a_mod_p);

        SetCoeff(smoothPair_mod_f_p._zz_pE__rep, 0, cd_x_a);
        SetCoeff(smoothPair_mod_f_p._zz_pE__rep, 1, b_mod_p);


        //conv(smoothPair_mod_f_p, smoothPair_mod_p);
        mul(dp, dp, smoothPair_mod_f_p);
    }


    zz_pX diff_g_mod_p;
    diff(diff_g_mod_p, g_mod_p);


    zz_pE diff_g_mod_f_p, diff_g_mod_f_p_squared;
    conv(diff_g_mod_f_p, diff_g_mod_p);

    mul(diff_g_mod_f_p_squared, diff_g_mod_f_p, diff_g_mod_f_p);


    mul(dp, dp, diff_g_mod_f_p_squared);


#ifdef GNFS_DEBUG
    assert(IsOne(power(dp, q_minus_1 >> 1)));
#endif

    return true;
}

inline void SquareRoot::calculate_r_s(long& r, ZZ& s) { //find r,s such that q - 1 = s * 2^r

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "==> In calculate_r_s()" << endl;
#endif
    s = q_minus_1;
    r = MakeOdd(s);

#ifdef GNFS_SQRT_STEP_DEBUG
    assert(q_minus_1 == s * power(to_ZZ(2), r));
#endif
}

bool SquareRoot::findQuadraticNonResidue(zz_pE& QNR) {
#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "==> In findQuadraticNonResidue()" << endl;
#endif

    static ZZ exp;

    RightShift(exp, q_minus_1, 1); //exp = q_minus_1 /2;

    zz_pE result;

    zz_pX poly_mod_p;
    SetCoeff(poly_mod_p, 1, 1);

    zz_p expectedResult;
    conv(expectedResult, -1);

    for (int i = 1; i < 100; i++) {

        SetCoeff(poly_mod_p, 0, i);
        //SetCoeff(poly_mod_p, 1, 1);

        conv(QNR, poly_mod_p);
        //cout << "polyext4 -> " << polyext4 << endl;

        power(result, QNR, exp);

        //cout << "dokimi -> " << dokimi << endl;
        if (deg(result._zz_pE__rep) == 0 && result._zz_pE__rep.rep[0] == expectedResult) {
            //cout << "QNR -> " << QNR << endl;

            return true;
        }
    }
    cout << "WARNING!: Could not find a quadratic non residue while processing prime " << prime << endl;
    return false;
}

/*
 * To calculate the square roots in the Finite Field defined from the selected prime,
 * we use the iterative approach described by Shanks & Tonelli. This approach is named
 * Shanks - Tonelli Algorithm.
 */
bool SquareRoot::calculateSquareRootOfFiniteField(zz_pE& root, const zz_pE& dp) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "==> In calculateSquareRootOfFiniteField()" << endl;
#endif

    /*
     * NOTE!!
     *
     * Here "root" is "w" as described in Shanks and Tonelli algorithm
     */

    static ZZ tmpZZ, tmpPow, s;
    long r;

    //calculate the values of r and s
    calculate_r_s(r, s);

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "r = " << r << endl;
    cout << "s = " << s << endl;
#endif

    zz_pE l; //l0,l1,l2 etc

    power(l, dp, s);
    //cout << "l0 -> " << l << endl;

    add(tmpZZ, s, 1);
    tmpZZ >>= 1; //tmpZZ == (s + 1) / 2
    power(root, dp, tmpZZ);



    //Find a Quadratic Non - Residue (QNR)
    zz_pE QNR;
    if (!findQuadraticNonResidue(QNR))
        return false;

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "QNR -> " << QNR << endl;
#endif

    zz_pE z, tmp_z;
    power(z, QNR, s);


    long m = 0;
    zz_pE power_of_l;

    while (true) {

        //We calculate the new m value
        //long order = 1;
        power_of_l = l;
        for (m = 0; m < r; ++m) {
            // power(power_of_l, l, order);
            if (IsOne(power_of_l))
                break;
            //order *= 2;
            mul(power_of_l, power_of_l, power_of_l);
        }


        if (m == 0) //we found the root! break!
            break;


        //If we did not find a proper m something must have gone wrong, so we return false...
        if (m == r) {
            cout << "ERROR!: Could not find a proper m while processing prime " << prime << endl;
            return false;
        }


        //compute the square root w(i) or "root" (from parameters)
        //w(i+1) = w(i)*z^(2^(r-m-1))
        power2(tmpPow, r - m - 1);
        power(tmp_z, z, tmpPow);
        mul(root, root, tmp_z);


        //compute l(i)
        //l(i+1) = l(i)*z^(2^(r-m))
        mul(z, tmp_z, tmp_z); //update z as well
        mul(l, l, z); //


        //If l == 1 we 're done
        if (IsOne(l)) {
            //if we did not break here, m would be found equal to 0 during the next iteration which is also a terminating condition
            break;
        }

        r = m;
    }

#ifdef GNFS_SQRT_STEP_DEBUG
    assert(root * root == dp);
#endif

    return true;
}

inline void SquareRoot::calculateNormBp(zz_p& norm_mod_p, const zz_pE& squareRoot) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "==> In calculateNormBp()" << endl;
#endif
    static ZZ exp;

    div(exp, q_minus_1, prime - 1);
    zz_pE norm_mod_f_p;
    power(norm_mod_f_p, squareRoot, exp);

    norm_mod_p = norm_mod_f_p._zz_pE__rep.rep[0];

#ifdef GNFS_DEBUG
    assert(deg(norm_mod_f_p._zz_pE__rep) == 0);
#endif
}

void SquareRoot::calculateNormBModPrime(zz_p& norm_mod_p, const zz_pX& f_mod_p, const long numOfElementsInDependency) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "==> In calculateNormBModPrime()" << endl;
#endif
    zz_p c_d_pow;
    power(c_d_pow, f_mod_p.rep[d], (d - 1) * (numOfElementsInDependency / 2));

    conv(norm_mod_p, squareRootOfNormOfDep_A);
    mul(norm_mod_p, norm_mod_p, c_d_pow);
    mul(norm_mod_p, norm_mod_p, to_zz_p(discr_g));

    // normN_a = to_zz_pE(norm);
    //cout << "N(a) -> " << normN_a << endl;

}

inline void SquareRoot::evaluate_xi(zz_p& xi, const zz_pE& squareRoot, const zz_pX& f_mod_p) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "==> In evaluate_xi()" << endl;
#endif
    zz_p x;
    conv(x, m);
    mul(x, x, f_mod_p.rep[d]);

    eval(xi, squareRoot._zz_pE__rep, x);
    //cout << "xi -> " << xi << endl;
}

void SquareRoot::crt(ZZ_p& x_mod_n, const ZZ_p& P_mod_n) {



    //We calculate "r"as described in Couveignes algorithm
    //  r : z/M = Sum( (ai*xi) / pi )
    ZZ_p r_mod_n;
    calculate_r_mod_n(r_mod_n);

    //
    //  X mod(n) = Sum( ai*xi*Pi( mod(n) ) ) - r*P( mod(n) )
    //
    calculate_x_mod_n(x_mod_n, r_mod_n, P_mod_n);
    //cout << "resXmodN -> " << resXmodN << endl;


}

/*void SquareRoot::calculate_P(ZZ& P) { //useless

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "=> In calculateP()" << endl;
#endif

    set(P);

    for (vector<long>::const_iterator iter = IPB.begin(); iter != IPB.end(); ++iter) {
        mul(P, P, *iter);
    }
}*/

void SquareRoot::calculate_P_mod_n(ZZ_p& P_mod_n) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "=> In calculate_P_mod_n()" << endl;
#endif

    set(P_mod_n);

    for (vector<long>::const_iterator iter = IPB.begin(); iter != IPB.end(); ++iter) {
        mul(P_mod_n, P_mod_n, *iter);
    }
}

void SquareRoot::calculateVectorOf_Pi_mod_p(vector<long>& Pi_mod_p_vec) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "=> In calculateVectorOf_Pi_mod_p()" << endl;
#endif
    Pi_mod_p_vec.clear();
    Pi_mod_p_vec.reserve(IPBSize);

    for (vector<long>::const_iterator iter = IPB.begin(); iter != IPB.end(); ++iter) {

        Pi_mod_p_vec.push_back(calculate_Pi_mod_p(*iter));
    }
}

inline long SquareRoot::calculate_Pi_mod_p(const long prime) {

    double pinv = 1 / ((double) prime);
    long Pi = 1;
    for (vector<long>::const_iterator iter = IPB.begin(); iter != IPB.end(); ++iter) {
        if (*iter != prime)
            Pi = MulMod(Pi, *iter, prime, pinv);
    }

    return Pi;
}

void SquareRoot::calculateVectorOf_ai(const vector<long>& Pi_mod_p_vec) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "=> In calculateVectorOf_ai()" << endl;
#endif

    ai_vec.clear();
    ai_vec.reserve(IPBSize);
    long i = 0;
    for (vector<long>::const_iterator iter = IPB.begin(); iter != IPB.end(); ++iter) {

        ai_vec.push_back(calculate_ai(Pi_mod_p_vec[i++], *iter));
    }
}

inline long SquareRoot::calculate_ai(const long Pi, const long prime) {

    return InvMod(Pi, prime);
}

void SquareRoot::calculate_r_mod_n(ZZ_p& r_mod_n) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "=> In calculate_r_mod_n()" << endl;
#endif
    long oldPrec = RR::precision();
    RR::SetPrecision(256L);

    ZZ r;
    RR sum;
    RR tmpTerm;

    long i = 0;
    for (vector<long>::const_iterator iter = IPB.begin(); iter != IPB.end(); ++iter) {

        /*conv(tmpProd,xi_vec[i]);
        mul(tmpProd,tmpProd,ai_vec[i]);*/
        conv(tmpTerm, xi_vec[i]);
        mul(tmpTerm, tmpTerm, ai_vec[i]);
        div(tmpTerm, tmpTerm, (*iter));
        // tmpTerm /= (*iter);
        add(sum, sum, tmpTerm);

        i++;
    }

    r = FloorToZZ(0.5 + sum);
    conv(r_mod_n, r);

    RR::SetPrecision(oldPrec);
    // cout << _r << endl;
    //  cout << r_mod_n << endl;
}

void SquareRoot::calculate_x_mod_n(ZZ_p& x_mod_n, const ZZ_p& r_mod_n, const ZZ_p & P_mod_n) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "=> In calculate_x_mod_n()" << endl;
#endif

    ZZ_p sumOfProduct_mod_n;
    ZZ_p r_x_P_mod_n;

    calculateSumOfProduct_ai_xi_Pi_mod_n(sumOfProduct_mod_n, P_mod_n);
    calculateProductOf_r_P_mod_n(r_x_P_mod_n, r_mod_n, P_mod_n);

    sub(x_mod_n, sumOfProduct_mod_n, r_x_P_mod_n);
}

void SquareRoot::calculateSumOfProduct_ai_xi_Pi_mod_n(ZZ_p& sumOfProduct_mod_n, const ZZ_p& P_mod_n) {

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "==> In calculateSumOfProduct_ai_xi_Pi_mod_n()" << endl;
#endif

    clear(sumOfProduct_mod_n);

    ZZ_p Pi_mod_n, tmpProd;

    long i = 0;
    for (vector<long>::const_iterator iter = IPB.begin(); iter != IPB.end(); ++iter) {

        div(Pi_mod_n, P_mod_n, *iter);

        mul(tmpProd, Pi_mod_n, ai_vec[i]);
        mul(tmpProd, tmpProd, xi_vec[i++]);

        add(sumOfProduct_mod_n, sumOfProduct_mod_n, tmpProd);
    }


}

inline void SquareRoot::calculateProductOf_r_P_mod_n(ZZ_p& r_x_P_mod_n, const ZZ_p& r_mod_n, const ZZ_p& P_mod_n) {

    mul(r_x_P_mod_n, r_mod_n, P_mod_n);
}

void SquareRoot::findSquareRootInZ_mod_n(ZZ_p& squareRoot_mod_n, const vec_Relation & dependency) {

    ZZ_p squareRootOfDep_R_mod_n;
    calculateSquareRootOfDepInZ_mod_n(squareRootOfDep_R_mod_n, dependency);

    ZZ_pX f_mod_n;
    conv(f_mod_n, f);

    ZZ_p derivSqrt_mod_n;
    calculate_diff_f_mod_n_in_m(derivSqrt_mod_n, f_mod_n);
    //cout << "derivSqrt -> " << derivSqrt << endl;

    long numOfElementsInDep = dependency.length();

    ZZ_p a_d_pow_mod_n;
    power(a_d_pow_mod_n, f_mod_n.rep[d], d - 2 + (numOfElementsInDep / 2));

    mul(squareRoot_mod_n, derivSqrt_mod_n, a_d_pow_mod_n);
    mul(squareRoot_mod_n, squareRoot_mod_n, squareRootOfDep_R_mod_n);

}

void SquareRoot::calculateSquareRootOfDepInZ_mod_n(ZZ_p& squareRoot_mod_n, const vec_Relation & dependency) {

    set(squareRoot_mod_n);

#ifdef GNFS_SQRT_STEP_DEBUG
    ZZ prod;
    set(prod);
#endif

    bzero(factor_exps_R, RFBSize * sizeof (unsigned short));
    for (long i = 0; i < dependency.length(); ++i) {
        const Relation& rel = dependency[i];

#ifdef GNFS_SQRT_STEP_DEBUG
        long a = dependency[i].relPair.a;
        long b = dependency[i].relPair.b;
        prod *= a - (b * m);
#endif

        long factorsLength = rel.rationalFactors.length();
        for (long k = 0; k < factorsLength; ++k) {
            factor_exps_R[rel.rationalFactorsIndex[k]] += rel.rationalFactors[k].b;
        }
    }

    ZZ tmpCoFactor;
    ZZ_p tmpCoFactor_mod_n;

    for (long k = 0; k < RFBSize; ++k) {
        power(tmpCoFactor, RFB[k].a, factor_exps_R[k] / 2);
        conv(tmpCoFactor_mod_n, tmpCoFactor);
        mul(squareRoot_mod_n, squareRoot_mod_n, tmpCoFactor_mod_n);
    }


#ifdef GNFS_SQRT_STEP_DEBUG
    SqrRoot(prod, prod);
    assert(prod % n == squareRoot_mod_n._ZZ_p__rep);
#endif
}

inline void SquareRoot::calculate_diff_f_mod_n_in_m(ZZ_p& derivSqrt_mod_n, const ZZ_pX & f_mod_n) {

    ZZ_pX diff_f_mod_n;
    diff(diff_f_mod_n, f_mod_n);
    eval(derivSqrt_mod_n, diff_f_mod_n, to_ZZ_p(m));

#ifdef GNFS_SQRT_STEP_DEBUG
    cout << "derivSqrt -> " << derivSqrt << endl;
#endif
}

void SquareRoot::calculate_gx(ZZX & g) { //calculates g(x) = f(x/c_d)*(c_d^(d - 1)), where c_d is the leading coefficient of f

    long d_minus_1 = d - 1;

    ZZ newCofactor;
    ZZ coefficient;
    long exp;

    g.kill();
    for (long i = 0; i < d; i++) {

        exp = d_minus_1 - i;
        power(newCofactor, f.rep[d], exp);

        mul(coefficient, newCofactor, f.rep[i]);

        SetCoeff(g, i, coefficient);
    }

    SetCoeff(g, d, 1L);
}

long SquareRoot::checkResult(ofstream& ofstr, const long currentDep, const ZZ& factor1, const ZZ& factor2) {

    if (IsOne(factor1) && IsOne(factor2)) {
        outputFailureMsg1(cout);
        outputFailureMsg1(ofstr);
        return -1;
    }

    cout << "\n*** Factorization result of " << n << endl << endl;
    cout << "The 1st factor is -> " << factor1 << endl;
    cout << "The 2nd factor is -> " << factor2 << endl;

    ofstr << "\n*** Factorization result of " << n << endl << endl;
    ofstr << "The 1st factor is -> " << factor1 << endl;
    ofstr << "The 2nd factor is -> " << factor2 << endl;

    ZZ check;
    mul(check, factor1, factor2);

    if (check == n) {

        if (IsOne(factor1) || IsOne(factor2)) {
            cout << endl << "Oops... Trivial factorization of n found..." << endl;
            ofstr << endl << "Oops... Trivial factorization of n found..." << endl;

            if (currentDep < dependencies.length() - 1) {

                cout << "The process continues to the next dependency..." << endl;
                ofstr << "The process continues to the next dependency..." << endl;
            }
            else {

                cout << "No more dependencies to process! Try increasing the number of fetched dependencies..." << endl;
                ofstr << "No more dependencies to process! Try increasing the number of fetched dependencies..." << endl;
                return -1;
            }
            cout << "======================================================" << endl;
            ofstr << "======================================================" << endl;
            return 0;
        }

        cout << "\n### SUCCESS!! ###" << endl;
        if (currentDep < 0)
            cout << endl << "Result found during the inert prime search" << endl << endl;
        else
            cout << endl << "Result found using dependency " << currentDep + 1 << " of " << dependencies.length() << endl << endl;
        cout << "*** Factorization of " << n << " completed successfully! " << endl;
        cout << "======================================================" << endl;

        ofstr << "\n### SUCCESS!! ###" << endl;
        if (currentDep < 0)
            ofstr << endl << "Result found during the inert prime search" << endl << endl;
        else
            ofstr << endl << "Result found using dependency " << currentDep + 1 << " of " << dependencies.length() << endl << endl;
        ofstr << "*** Factorization of " << n << " completed successfully! " << endl;
        ofstr << "======================================================" << endl;

        return 1;
    }

    outputFailureMsg2(cout);
    outputFailureMsg2(ofstr);
    return -1;
}

void SquareRoot::outputFailureMsg1(ostream& ostr) {

    ostr << "### FAILURE ###" << endl;
    ostr << "Couveignes square root step failed to work." << endl;
    ostr << "Both GCDs returned 1 thus meaning that more inert primes pi are necessary in order to make a successful calculation." << endl;
    ostr << "Try increasing the number of inert primes and/or their size through initialization_parameters.txt" << endl;
}

void SquareRoot::outputFailureMsg2(ostream& ostr) {

    ostr << "### FAILURE ###" << endl;
    ostr << "The result is not correct...! This should not happen........" << endl;
}

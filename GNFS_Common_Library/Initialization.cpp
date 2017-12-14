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


#include <GNFS/Initialization.h>


const char* Initialization::inputParamsFile = "IO/Initialization/initialization_parameters.txt"; //The initial parameters file

void Initialization::initialCalculations() {

    if (NumBits(n) > RR::prec - 10)
        precision = NumBits(n) + 10;
    else
        precision = RR::prec;

    if (m == 0) {
        calculateMinM();
        calculateMaxM();
        m = maxM;
    }
    else {
        maxM = minM = m;
    }

    if (calculateRFBAndAFBBound) {
        cout << "Programmatically calculating factor base bounds..." << endl;
        calculateBounds();
    }

    if (L1 <= maxRFB || L2 <= maxAFB) {
        L1 = maxRFB;
        L2 = maxAFB;
    }


    f = ZZX::zero();
    diffF = ZZX::zero();

}

void Initialization::calculateMinM() {

    RR::SetPrecision(precision);
    FloorToZZ(minM, pow(to_RR(n), to_RR(1.0) / (d + 1)));


}

void Initialization::calculateMaxM() {

    RR::SetPrecision(precision);
    FloorToZZ(maxM, pow(to_RR(n), to_RR(1.0) / d));

}

void Initialization::calculateBounds() {//Better not to use this. Enter the bounds manually, based on previous experience.

    //This formula is from Pomerance's book "Prime Numbers: A Computational Perspective", chapter six.
    //It seems to yield much higher bounds than those that are used in practice
    RR::SetPrecision(precision);
    RR n = to_RR(this->n);
    RR int1 = pow(to_RR(8.0) / 9.0, to_RR(1.0) / 3.0);
    RR int2 = pow(log(n), to_RR(1.0) / 3.0);
    RR int3 = pow(log(log(n)), to_RR(2.0) / 3.0);

    conv(maxRFB, FloorToZZ(exp(int1 * int2 * int3)));
    maxAFB = maxAFBFactor*maxRFB;

    if (L1 != 0)
        L1 = static_cast<long> (floor(pow(maxRFB, LExponent))); //todo future int64_type
    if (L2 != 0)
        L2 = static_cast<long> (floor(pow(maxAFB, LExponent)));

}

void Initialization::calculateD() {//Not in use, d is read from the input file directly

    RR::SetPrecision(precision);
    RR n = to_RR(this->n);

    RR med = 3 * log(n) / log(log(n));

    d = to_long(FloorToZZ(pow(med, to_RR(1.0) / 3.0)));
}

bool Initialization::runInitialTests() const {

    cout << "### Please read carefully the following system parameters!!! ###" << endl;
    cout << "-> char: " << sizeof (char) << " bytes..." << endl;
    cout << "-> short: " << sizeof (short) << " bytes..." << endl;
    cout << "-> int: " << sizeof (int) << " bytes..." << endl;
    cout << "-> long: " << sizeof (long) << " bytes..." << endl;
    cout << "-> int_type: " << sizeof (int_type) << " bytes..." << endl;
    cout << "-> uint_type: " << sizeof (uint_type) << " bytes..." << endl;
    cout << "-> NTL single precision bound (NTL_SP_BOUND): \t\t" << NTL_SP_BOUND << endl;
    cout << "-> NTL single precision bits (NTL_SP_NBITS): " << NTL_SP_NBITS << " bits..." << endl;
    cout << "-> NTL wide single precision bound (NTL_WSP_BOUND): \t" << NTL_WSP_BOUND << endl;
    cout << "-> NTL wide single precision bits (NTL_SWP_NBITS): " << NTL_WSP_NBITS << " bits..." << endl;
    cout << "-> NTL bits per zzigit (NTL_ZZ_NBITS): " << NTL_ZZ_NBITS << " bits..." << endl;
    cout << "-> NTL bits per int (NTL_BITS_PER_INT): " << NTL_BITS_PER_INT << " bits..." << endl;
    cout << "-> NTL bits per long (NTL_BITS_PER_LONG): " << NTL_BITS_PER_LONG << " bits..." << endl;

    if (NTL_BITS_PER_LONG == 64)
        cout << "##### Good news! Looks like we have a 64-bit system! #####" << endl;
    else
        cout << "@@@@@ Warning! Looks like we have a 32-bit system! long is only 4 bytes. Be careful... @@@@@" << endl;

    if (sizeof (int_type) != sizeof (uint_type)) {
        cout << "int_type and uint_type must have the same size!" << endl;
        return false;
    }
    if (sizeof (int_type) > sizeof (long)) {
        cout << "The size of int_type and uint_type must no be larger than the size of long!" << endl;
        return false;
    }

    /*if (sizeof (int64_type) != 8 || sizeof (uint64_type) != 8) {
        cout << "Please define trully 64-bit data types as int64_type and uint64_type in global.h, depending on the current platrform, then rebuild and run again! Aborting..." << endl;
        return false;
    }*///todo future

    if (maxAFB % 2 != 0 || maxRFB % 2 != 0 || L1 % 2 != 0 || L2 % 2 != 0) {
        cout << "Please select maxAFB, maxRFB, L1 and L2 as even numbers!" << endl;
        return false;
    }

    if (isNPrime()) {
        cout << "Primality test for n showed that n is prime..." << endl;
        return false;
    }

    if (!testD()) {
        cout << "Test for d failed..." << endl;
        return false;
    }

    testM();

    return true;

}

bool Initialization::isNPrime() const {//Check that n is not prime. This is the only check performed on n. No check is performed to determine whether n is power of a prime p.

    return (ProbPrime(n, numberOfMillerRabinTests) == 0) ? false : true;

}

bool Initialization::testD() const {//Check that d^(2d^2) < n

    long exp = 2L * power_long(d, 2L);
    ZZ result = power_ZZ(d, exp);
    //cout << result << endl;
    if (result < n)
        return true;

    if (d == 3) {
        cout << "Warning: Test for d failed, but d == 3, which is the smallest acceptable value, so the test is considered successful..." << endl;
        return true;
    }
    return false;

}

bool Initialization::testM() const {//Check that gcd(m,n) = 1. If not, which is higly improbable, then we have factored n!

    ZZ g = GCD(n, m);
    if (g != 1) {
        cout << "Impossible luck! n and m share a common factor: " << g << ". n is factored!" << endl;
        exit(EXIT_SUCCESS);
    }

    return true;
}

bool Initialization::polynomialSelectionStep() {
    PolynomialSelection polySel(n, d);
    if (monic) {

        if (extensivePolynomialSelection) {
            cout << "Finding the best rated monic polynomial through an extensive search..." << endl;

            if (byRating)
                cout << "*** Polynomial selection is performed by the highest rating value..." << endl;
            else
                cout << "*** Polynomial selection is performed by the least alpha value..." << endl;

            bool performSieveTest = false; //todo --- The "true" case is to be implemented, and we should read it from parameter file
            long numOfBestRatedPolynomialsToStore = 50; //todo ... read from parameter file
            //bool ratingWithPrecision_RR = false; //todo ... read from parameter file


            /*
             * NOTE:
             *
             * The polynomial selection with the RR precision is commented and never used
             * in practise, because the results are adequate with just double precision,
             * which is also faster. The commented code is not deleted, just in case we
             * want to test sth in the near future.
             */

            /*
                  if (ratingWithPrecision_RR) {
                    bool res = polySel.findBestPolynomial_RR(f, m, minM, maxM, maxAFB, numOfBestRatedPolynomialsToStore, performSieveTest);
                    if (!res)
                        return false;
                    }
                  else {
             */

            //bool res = polySel.findBestPolynomial(f, m, minM, maxM, maxAFB, numOfBestRatedPolynomialsToStore, performSieveTest);
            bool res = polySel.findBestPolynomial(f, m, minM, maxM, maxAFB, byRating, numOfBestRatedPolynomialsToStore, performSieveTest);
            if (!res)
                return false;

            /*
                } //EOF if (ratingWithPrecision_RR) {}
             */
        }
        else {
            cout << "Finding a monic polynomial without any optimization..." << endl;
            f = polySel.findMonicPolynomial(m);
        }
    }
    else {
        if (extensivePolynomialSelection) {
            cout << "Finding non-monic, non-skewed polynomial through an extensive search..." << endl;

            if (byRating)
                cout << "*** Polynomial selection is performed by the highest rating value..." << endl;
            else
                cout << "*** Polynomial selection is performed by the least alpha value..." << endl;

            //todo
            bool performSieveTest = false; //todo --- The "true" case is to be implemented, and we should read it from parameter file
            //ZZ cofactor = to_ZZ(2 * 2 * 3 * 5); //todo --- read from parameter file. It should be product of small primes, i.e 2*3*5*7*11
            long numOfBestRatedPolynomialsToStore = 20; //todo ... read from parameter file

            polySel.findBestNonSkewedNonMonicPolynomial(f, m, maxAFB, minM, maxM, x1, x2, range, byRating, numOfBestRatedPolynomialsToStore, performSieveTest);
        }
        else {
            cout << "Finding a non-monic, non-skewed polynomial without any optimization..." << endl;

            double logN = log(n);
            double x1 = 0.02;
            ZZ leadingCoef = to_ZZ(exp((d * log(x1) + logN) / (d + 1)));


            //todo
            //
            //Na ftiaksw to precision otan kalw RR
            RR tempPow;
            pow(tempPow, to_RR(n) / to_RR(leadingCoef), to_RR(1) / d);
            FloorToZZ(m, tempPow);

            //cout << endl << "m selected -> " << m << endl << endl;

            f = polySel.findNonSkewedNonMonicPolynomialWithLeadingCoefKnown(m, leadingCoef);

            vector<long> primesLessThan100;
            polySel.getPrimesLessThan100(primesLessThan100);

            double alpha = polySel.calcAlpha(f, primesLessThan100);
            double rating = polySel.ratePolynomial(f, alpha, maxAFB);

            cout << "---------------------------------------------------------" << endl;
            cout << "---Non-skewed, Non-monic polynomial selection results ---" << endl;
            cout << "---------------------------------------------------------" << endl;
            cout << "-- polynomial -> " << f << endl;
            cout << "-- m -> " << m << endl;
            cout << "-- alpha -> " << alpha << endl;
            cout << "-- rating -> " << rating << endl;
            cout << "---------------------------------------------------------" << endl;
        }
    }

    cout << "Polynomial selection finished..." << endl << "The selected polynomial is f = " << f << ", m = " << m << endl << "Trying to factor f in Z[x]...";
    if (verbose_ZZX_Factoring)
        cout << "Verbose mode...";

    cout << endl;
    vec_pair_ZZX_long factors;
    ZZ c;

    ZZXFac_InitNumPrimes = InitNumPrimes_ZZX_Factoring;
    factor(c, factors, f, verbose_ZZX_Factoring);


    if (factors.length() == 1 && c == 1)
        cout << "f is irreducible..." << endl;
    else if (factors.length() == 1 && c != 1) {
        cout << "Impossible luck! f is reducible! This yields a non trivial factorization of n!" << endl;
        cout << "Printing out the factorization result..." << endl;
        cout << "Content of f: " << c << endl;
        cout << "Factors of f: " << factors << endl;
        printFactorizationOfN(n, c);
        return true;
    }
    else if (factors.length() == 2 && c == 1) {
        cout << "Impossible luck! f is reducible! This yields a non trivial factorization of n!" << endl;
        cout << "Printing out the factorization result..." << endl;
        cout << "Content of f: " << c << endl;
        cout << "Factors of f: " << factors << endl;
        ZZX g = factors[0].a;
        ZZ factorOfN;
        eval(factorOfN, g, m);
        printFactorizationOfN(n, factorOfN);
        return true;
    }
    else {
        cout << "Impossible luck! f has more than 2 factors in Z[x]. This shows that n has more than 2 factors. Printing out the factorization result..." << endl;
        cout << "Content of f: " << c << endl;
        cout << "Factors of f: " << factors << endl;
        return true;
    }

    cout << "Calculating and testing the derivative of f..." << endl;
    diffF = diff(f);

    ZZ g = GCD(eval(diffF, m), n);
    if (g != 1 && g != n) {
        cout << "GCD(f'(m), n) is a non trivial factor of n!" << endl;
        printFactorizationOfN(n, g);
        return true;
    }

    GCD(g, coeff(f, d), n);
    if (g != 1 && g != n) {
        cout << "GCD(c_d, n) is a non trivial factor of n!" << endl;
        printFactorizationOfN(n, g);
        return true;
    }


    return true;
}

void Initialization::readInputFileStream(ifstream& f) {

    Parameters::readInputFileStream(f);

    f >> boolalpha >> performPolynomialSelection;
    f.ignore(INT_MAX, '\n');
    f >> boolalpha >> extensivePolynomialSelection;
    f.ignore(INT_MAX, '\n');

    //todo
    f >> x1;
    f.ignore(INT_MAX, '\n');
    f >> x2;
    f.ignore(INT_MAX, '\n');
    f >> range;
    f.ignore(INT_MAX, '\n');

    f >> boolalpha >> byRating;
    f.ignore(INT_MAX, '\n');

    f >> boolalpha >> calculateRFBAndAFBBound;
    f.ignore(INT_MAX, '\n');

    if (x2 < x1) {
        cout << "Warning: x1 should be smaller than x2... Fixing this by swapping the two values..." << endl;
        //x1 = x2;
        swap(x1, x2);
    }
}

void Initialization::writeOutputFileStream(ofstream& f) const {

    Parameters::writeOutputFileStream(f);

    f << "# ADDITIONAL NOTES:" << endl;
    if (!performPolynomialSelection)
        f << "# The polynomial selection step was skipped. No polynomial f was selected." << endl;
    f << "# ATTENTION!!: Large prime optimization is not implemented yet!!" << endl; //todo future

}
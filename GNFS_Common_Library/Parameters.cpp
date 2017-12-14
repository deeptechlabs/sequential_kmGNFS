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


#include <GNFS/Parameters.h>

const char* Parameters::paramsFile = "IO/Initialization/gnfs_parameters.txt"; //gnfs algorithm parameters

void Parameters::readInputFileStream(ifstream& istr) {

    istr >> n;
    istr.ignore(INT_MAX, '\n');
    istr >> d;
    istr.ignore(INT_MAX, '\n');
    istr >> boolalpha >> monic;
    istr.ignore(INT_MAX, '\n');
    istr >> f;
    istr.ignore(INT_MAX, '\n');
    istr >> diffF;
    istr.ignore(INT_MAX, '\n');
    istr >> m;
    istr.ignore(INT_MAX, '\n');
    istr >> maxRFB;
    istr.ignore(INT_MAX, '\n');
    istr >> maxAFB;
    istr.ignore(INT_MAX, '\n');
    istr >> QCBSize;
    istr.ignore(INT_MAX, '\n');
    istr >> boolalpha >> outputBasesInText;
    istr.ignore(INT_MAX, '\n');
    istr >> max_a;
    istr.ignore(INT_MAX, '\n');
    istr >> min_b;
    istr.ignore(INT_MAX, '\n');
    istr >> max_b;
    istr.ignore(INT_MAX, '\n');
    istr >> L1;
    istr.ignore(INT_MAX, '\n');
    istr >> L2;
    istr.ignore(INT_MAX, '\n');
    istr >> numOfExcessiveRelations;
    istr.ignore(INT_MAX, '\n');
    istr >> minInertPrimeLength;
    istr.ignore(INT_MAX, '\n');
    istr >> IPBSize;
    istr.ignore(INT_MAX, '\n');
    istr >> maxNumOfDepsToProcess;
    istr.ignore(INT_MAX, '\n');

    if (maxAFB < maxRFB) {
        cout << "Warning: maxRFB should not be larger than maxAFB... Fixing this by setting maxRFB = maxAFB..." << endl;
        maxRFB = maxAFB;
    }

    if (L1 <= maxRFB || L2 <= maxAFB) {
        L1 = maxRFB;
        L2 = maxAFB;
    }

    if (L2 < L1) {
        cout << "Warning: L1 should not be larger than L2... Fixing this by setting L1 = L2..." << endl;
        L1 = L2;
    }

    if (min_b <= 0)
        min_b = 1;

    if (max_b < min_b)
        max_b = min_b;

    if (coeff(f, d) > 1) {//make sure monic flag is false if f is indeed non monic
        monic = false;
    }

    if (minInertPrimeLength > NTL_SP_NBITS) {
        cout << "Warning: The length in bits of the inert primes should not be greater than NTL_SP_NBITS (" << NTL_SP_NBITS << ") . Deacreasing its value to NTL_SP_NBITS..." << endl;
        minInertPrimeLength = NTL_SP_NBITS;
    }

    if (minInertPrimeLength < 14) { //fix minInertPrimeLength if it is too small
        cout << "Warning: The length in bits of the inert primes is too small. Fixing it to 14 bits..." << endl;
        minInertPrimeLength = 14;
    }

    if (IPBSize <= 0)
        IPBSize == 1;

    if (maxNumOfDepsToProcess <= 0)
        maxNumOfDepsToProcess = 4;

}

void Parameters::writeOutputFileStream(ofstream& ostr) const {

    ostr << n;
    ostr << "\t\t\t\t" << "# Number to factor (n)" << endl;
    ostr << d;
    ostr << "\t\t\t\t" << "# Polynomial degree (d)" << endl;
    ostr << boolalpha << monic;
    ostr << "\t\t\t\t" << "# If true, the polynomial that will be used by GNFS will be monic." << endl;
    ostr << f;
    ostr << "\t\t\t\t" << "# Selected polynomial (f). NOTE: This will be empty if the polynomial selection step was skipped." << endl;
    ostr << diffF;
    ostr << "\t\t\t\t" << "# Derivative of f. NOTE: This will be empty if the polynomial selection step was skipped." << endl;
    ostr << m;
    ostr << "\t\t\t\t" << "# Root of f mod n (m)" << endl;
    ostr << maxRFB;
    ostr << "\t\t\t\t" << "# Rational factor base bound (maxRFB)" << endl;
    ostr << maxAFB;
    ostr << "\t\t\t\t" << "# Algebraic factor base bound (maxAFB). It must be >= maxRFB" << endl;
    ostr << QCBSize;
    ostr << "\t\t\t\t" << "# Quadratic character base size (#QCB)" << endl;
    ostr << boolalpha << outputBasesInText;
    ostr << "\t\t\t\t" << "# false: Output bases only in binary format. true: Output bases in text format as well." << endl;
    ostr << max_a;
    ostr << "\t\t\t\t" << "# x dimension size (actually the half of it) of the line sieving step ( |a| <= max_a )" << endl;
    ostr << min_b;
    ostr << "\t\t\t\t" << "# lower bound for y dimension of the line sieving step (min_b). Must be at least one." << endl;
    ostr << max_b;
    ostr << "\t\t\t\t" << "# upper bound for y dimension of the line sieving step (max_b, min_b <= b <= max_b)" << endl;
    ostr << L1;
    ostr << "\t\t\t\t" << "# Large prime bound for RFB (L1). It must be lagrer than maxRFB, or else this optimization will not be applied." << endl;
    ostr << L2;
    ostr << "\t\t\t\t" << "# Large prime bound for AFB (L2). It must be lagrer than maxAFB, or else this optimization will not be applied. It must be >= L1" << endl;
    ostr << numOfExcessiveRelations;
    ostr << "\t\t\t\t" << "# Minimum number of excessive relations to incorporate into the matrix during the linear algebra step." << endl;
    ostr << minInertPrimeLength;
    ostr << "\t\t\t\t" << "# Starting length in bits of the inert primes that will be used in the square root step. The bigger this value is the fewer inert primes will be necessary (i.e. IPBsize can be smaller). See the comments in Parameters.h for more info!" << endl;
    ostr << IPBSize;
    ostr << "\t\t\t\t" << "# Number of inert primes that will be used for the square root step. See the comments in Parameters.h for more info!" << endl;
    ostr << maxNumOfDepsToProcess;
    ostr << "\t\t\t\t" << "# Maximum number of dependencies to process in the square root step. See the comments in Parameters.h for more info!" << endl;
}

void Parameters::readInputParamsFile(const char* inputFile) {

    ifstream f(inputFile);

    if (!f.is_open()) {
        cout << "Error opening " << inputFile << ". Exiting..." << endl;
        exit(EXIT_FAILURE);
    }
    readInputFileStream(f);
    f.close();
}

void Parameters::writeParamsFile() const {

    ofstream f(paramsFile);

    if (!f.is_open()) {
        cout << "Error opening " << paramsFile << ". Exiting..." << endl;
        exit(EXIT_FAILURE);
    }
    writeOutputFileStream(f);
    f.close();
}

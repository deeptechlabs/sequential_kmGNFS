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
#include <cmath>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/RR.h>

#include <GNFS/gnfs_utils.h>

NTL_CLIENT

        const long precision = 256; //for the floating point versions of evalHomogeneousF

void eval(ZZ& res, const ZZX& f, const ZZ& a) {//Horner evaluation for ZZX class, procedural version
    //CAUTION: res must not be an alias of a!!

    static ZZ tmp;

    clear(res);
    for (long i = deg(f); i >= 0; --i) {
        mul(tmp, res, a);
        add(res, tmp, f.rep[i]);
    }
}

void evalHomogeneousF(ZZ& res, const ZZX& f, const ZZ& x, const ZZ& y) {//Evaluation of the homogeneous polynomial F(x,y) = y^d * f(x/y) = c_d*x^d + c_d_1*x^(d-1)*y + ... + c_0*y
    //CAUTION: res must not be an alias of x or y!!

    static ZZ tmp, x_pow;

    long d = deg(f);

    set(x_pow); //x_pow = 1

    res = f.rep[0];

    for (long i = 1; i <= d; ++i) {
        mul(res, res, y);
        mul(x_pow, x_pow, x);
        mul(tmp, x_pow, f.rep[i]);
        add(res, res, tmp);
    }
}

double evalHomogeneousF(const ZZX& f, const double x, const double y) {//Evaluation of the homogeneous polynomial F(x,y) = y^d * f(x/y) = c_d*x^d + c_d_1*x^(d-1)*y + ... + c_0*y
    //Floting point version using doubles for input/output and RRs for the intermediate calculations

    static RR tmp, x_pow, RR_x, RR_y, res;

    long oldPrec = RR::precision();

    RR::SetPrecision(precision);

    long d = deg(f);

    conv(RR_x, x);
    conv(RR_y, y);


    set(x_pow); //x_pow = 1

    conv(res, f.rep[0]);

    for (long i = 1; i <= d; ++i) {
        mul(res, res, RR_y);
        mul(x_pow, x_pow, RR_x);
        conv(tmp, f.rep[i]);
        mul(tmp, x_pow, tmp);
        add(res, res, tmp);
    }

    RR::SetPrecision(oldPrec);

    return to_double(res);

}

void evalHomogeneousF(RR& res, const ZZX& f, const RR& x, const RR& y) {//Evaluation of the homogeneous polynomial F(x,y) = y^d * f(x/y) = c_d*x^d + c_d_1*x^(d-1)*y + ... + c_0*y
    //Floting point version using only RRs
    //CAUTION: res must not be an alias of x or y!!

    static RR tmp, x_pow;

    long oldPrec = RR::precision();

    if (oldPrec < precision)
        RR::SetPrecision(precision);

    long d = deg(f);


    set(x_pow); //x_pow = 1

    conv(res, f.rep[0]);

    for (long i = 1; i <= d; ++i) {
        mul(res, res, y);
        mul(x_pow, x_pow, x);
        conv(tmp, f.rep[i]);
        mul(tmp, x_pow, tmp);
        add(res, res, tmp);
    }

    RR::SetPrecision(oldPrec);

}

void printFactorizationOfN(const ZZ& n, const ZZ& aFactorOfN) {

    ZZ otherFactor = n / aFactorOfN;
    cout << "Factorization of n = " << n << endl;
    cout << "n = " << aFactorOfN << " * " << otherFactor << endl;

}

void printGNFSUsage() {

    cout << "------ USAGE ------" << endl;
    cout << "All GNFS subprograms take as a single command line argument the path " <<
            "to the I/O directory, i.e. the directory where all the input or " <<
            "output files that are used by the subprograms are stored." << endl <<
            "The I/O directory which by default has the name \"GNFS_IO_Files\" " <<
            "contains a directory under the name \"IO\" which in " <<
            "turn contains a number of directories that describe the basic steps " <<
            "of the algorithm and contain the appropriate I/O files." << endl <<
            "NOTE: The supplied path becomes immediately the current working directory " <<
            "of the running process." << endl <<
            "NOTE: The supplied path could be absolute or relative." << endl <<
            "EXAMPLE: ./gnfsprogram \"/path/to/IODir/GNFS_IO_FILES\" " << endl <<
            "------------------" << endl;
}

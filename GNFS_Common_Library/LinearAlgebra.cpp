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
#include <NTL/vec_vec_GF2.h>
#include <NTL/mat_GF2.h>


#include <iostream>
#include <string>
#include <cstdlib>
#include <climits>
#include <dirent.h>
#include <cerrno>
#include <cstring>

#include <GNFS/Relation.h>
#include <GNFS/gnfs_utils.h>
#include <GNFS/gnfs_global.h>

#include <GNFS/LinearAlgebra.h>

#include <timer.h>


const char* LinearAlgebra::reducedMatrixTextFile = "IO/Matrix/reduced_matrix.txt";
const char* LinearAlgebra::reducedMatrixBinaryFile = "IO/Matrix/reduced_matrix.bin";

NTL_CLIENT

bool LinearAlgebra::init() {

    if (!MatrixConstruction::init())
        return false;

    gnfsMatrixRank = 0;
    linearAlgebraInitialized = true;
    return true;
}

bool LinearAlgebra::reduceMatrixUsingGauss(bool readMatrixFromDisk, bool writeReducedMatrixToDisk, bool binaryIO) {

    if (readMatrixFromDisk) {
        gnfsMatrix.kill();

        cout << "Reading constructed matrix from disk..." << endl;
        if (binaryIO) {
            if (!MatrixConstruction::readMatrixBinary(matrixBinaryFile))
                return false;
        }
        else {
            if (!MatrixConstruction::readMatrixText(matrixTextFile))
                return false;
        }
    }

    cout << "Performing gaussian elimination... " << endl;
    gnfsMatrixRank = gauss();
    cout << "Matrix rank: " << gnfsMatrixRank << endl;

    if (writeReducedMatrixToDisk) {
        cout << "Writing reduced matrix to disk..." << endl;
        if (binaryIO)
            return writeMatrixBinary(reducedMatrixBinaryFile);
        else
            return writeMatrixText(reducedMatrixTextFile);
    }

    return true;
}

long LinearAlgebra::gauss() {

    mat_GF2 &M = gnfsMatrix;

    const long n = M.NumRows();
    const long m = M.NumCols();

    freeVariables.SetLength(0);
    freeVariables.SetLength(m); //this resizing effectively zeroes out the bit vector

    const long wm = (m + NTL_BITS_PER_LONG - 1) / NTL_BITS_PER_LONG;

    long k, l;
    long i, j;
    long pos;

    l = 0;
    for (k = 0; k < m && l < n; ++k) {
        const long wk = k / NTL_BITS_PER_LONG;
        const long bk = k - wk*NTL_BITS_PER_LONG;
        const _ntl_ulong k_mask = 1UL << bk;


        pos = -1;
        for (i = l; i < n; ++i) {
            if (M[i].rep.elts()[wk] & k_mask) {
                pos = i;
                break;
            }
        }

        if (pos != -1) {
            if (l != pos)
                swap(M[pos], M[l]);

            _ntl_ulong *y = M[l].rep.elts();

            for (i = l + 1; i < n; ++i) {
                // M[i] = M[i] + M[l]*M[i,k]

                if (M[i].rep.elts()[wk] & k_mask) {
                    _ntl_ulong *x = M[i].rep.elts();

                    for (j = wk; j < wm; ++j)
                        x[j] ^= y[j];
                }
            }

            ++l;
        }
        else {
            freeVariables.rep.elts()[wk] |= k_mask;

        }
    }

    if (k < m) {
        GF2 one;
        set(one);
        for (; k < m; ++k)
            freeVariables.put(k, one);
    }

#ifdef GNFS_DEBUG
    assert(l + weight(freeVariables) == m);
#endif

    return l; // l is the rank of the matrix
}

bool LinearAlgebra::getSolutions(mat_GF2& solutions, const long numberOfSolutions) {

    assert(numberOfSolutions <= (gnfsMatrix.NumCols() - gnfsMatrixRank));

    solutions.kill();
    solutions.SetDims(numberOfSolutions, gnfsMatrix.NumCols());

    for (long i = 0; i < numberOfSolutions; ++i) {
        if (!getNextSolution(solutions[i]))
            return false;
    }

    return true;
}

bool LinearAlgebra::getDependenciesFromSolutions(vec_vec_Relation& dependencies, const mat_GF2& solutions) {

    if (!initializeRelationsDirectory())
        return false;


    long numOfSolutions = solutions.NumRows(); //each row is a single solution

    assert(numOfRequiredRelations == solutions.NumCols());

    dependencies.SetLength(numOfSolutions);

    for (long n = 0; n < numOfSolutions; ++n)
        dependencies[n].SetLength(0);

    Relation rel;
    for (long j = 0; j < numOfRequiredRelations; ++j) {
        if (!getNextRelationFromDisk(rel))
            return false;

        long a = rel.relPair.a;
        long b = rel.relPair.b;

        long wj = j / NTL_BITS_PER_LONG;
        long bj = j - wj*NTL_BITS_PER_LONG;
        _ntl_ulong j_mask = 1UL << bj;

        for (long n = 0; n < numOfSolutions; ++n) {
            if (solutions[n].rep.elts()[wj] & j_mask)
                append(dependencies[n], rel);
        }

    }



    return true;
}

bool LinearAlgebra::getNextSolution(vec_GF2& solution) { //solution must be of the proper length  and filled with zeroes

    static long lastSetFreeVarIndex = -1; //holds the index of the most recent free variable that was set to 1

    const long columns = gnfsMatrix.NumCols();

    if (solution.length() != columns) {
        cout << "Error while trying to get a solution. The supplied argument is not of proper size..." << endl;
        return false;
    }


    long i;
    for (i = lastSetFreeVarIndex + 1; i < columns; ++i) {
        if (IsOne(freeVariables.get(i))) {
            solution.put(i, 1L);
            lastSetFreeVarIndex = i;
            break;
        }

    }
    if (i >= columns)
        return false;

    return getSolution(solution);
}

bool LinearAlgebra::getSolution(vec_GF2& solution) {

    const long rows = gnfsMatrix.NumRows();
    const long columns = gnfsMatrix.NumCols();
    GF2 one;
    set(one);

    vec_GF2 computedVariables = freeVariables;
    //long wm = (columns + NTL_BITS_PER_LONG - 1) / NTL_BITS_PER_LONG;

    for (long i = rows - 1; i >= 0; --i) {

        long rowSum = 0;
        for (long j = columns - 1; j >= 0; --j) {
            /*if (IsOne(gnfsMatrix.get(i, j))) {
                if (IsOne(computedVariables.get(j)))
                    rowSum ^= rep(solution.get(j));
                else {
                    computedVariables.put(j, one);
                    solution.put(j, rowSum);
                    break;
                }
            }*/

            //find the bit that corresponds to column j and create the appropriate mask for manipulating this bit
            long wj = j / NTL_BITS_PER_LONG;
            long bj = j - wj*NTL_BITS_PER_LONG;
            _ntl_ulong j_mask = 1UL << bj;

            if (gnfsMatrix[i].rep.elts()[wj] & j_mask) {
                if (computedVariables.rep.elts()[wj] & j_mask)
                    rowSum ^= ((solution.rep.elts()[wj] >> bj) & 1);
                else {
                    computedVariables.rep.elts()[wj] |= j_mask;
                    if (rowSum == 0)
                        solution.rep.elts()[wj] &= ~j_mask;
                    else
                        solution.rep.elts()[wj] |= j_mask;
                    break;
                }
            }

        }

    }

    return true;
}

bool LinearAlgebra::writeMatrixText(const char* filename) {

    ofstream f;
    f.open(filename);
    if (!f.is_open()) {
        cout << "Error opening " << filename << "..." << endl;
        return false;
    }

    f << gnfsMatrix << endl << freeVariables << endl << gnfsMatrixRank << endl;

    if (f.good()) {
        f.close();
        return true;
    }

    cout << "Error while reading from " << filename << "..." << endl;
    f.close();
    return false;
}

bool LinearAlgebra::writeMatrixBinary(const char* filename) {

    if (!MatrixConstruction::writeMatrixBinary(filename))
        return false;

    long columns = freeVariables.length(); //this should be equal to gnfsMatrix.NumCols() or 0 (uninitialized)

    assert(columns == 0 || columns == gnfsMatrix.NumCols());

    ofstream f;
    f.open(filename, ios::out | ios::app | ios::binary);
    if (!f.is_open()) {
        cout << "Error opening " << filename << "..." << endl;
        return false;
    }

    if (columns != 0) {
        long wordsLength = (columns / NTL_BITS_PER_LONG) + 1;
        for (long j = 0; j < wordsLength; ++j)
            f.write(reinterpret_cast<char*> (freeVariables.rep.elts() + j), sizeof (_ntl_ulong)); //write freeVariables data
    }
    f.write(reinterpret_cast<char*> (&columns), sizeof (long)); //write freVariables columns number
    f.write(reinterpret_cast<char*> (&gnfsMatrixRank), sizeof (long)); //write rank of reduced matrix

    if (!f.good()) {
        f.close();
        return false;
    }

    f.close();
    return true;
}

bool LinearAlgebra::readMatrixText(const char* filename) {

    ifstream f;
    f.open(filename);
    if (!f.is_open()) {
        cout << "Error opening " << filename << "..." << endl;
        return false;
    }

    f >> gnfsMatrix >> freeVariables >> gnfsMatrixRank;

    if (f.good()) {
        f.close();
        return true;
    }

    cout << "Error while reading from " << filename << "..." << endl;
    f.close();
    return false;
}

bool LinearAlgebra::readMatrixBinary(const char* filename) {

    if (!MatrixConstruction::readMatrixBinary(filename))
        return false;

    ifstream f;
    f.open(filename, ios::in | ios::binary);
    if (!f.is_open()) {
        cout << "Error opening " << filename << "..." << endl;
        return false;
    }

    long columns = 0; //columns should be equal to gnfsMatrix.NumCols() or 0


    //Read one long before the last to figure out the length (columns) of freeVariables
    f.seekg(-2 * sizeof (long), ios::end);
    f.read(reinterpret_cast<char*> (&columns), sizeof (long));
    //Read the last long to figure out the rank of the reduced matrix
    f.read(reinterpret_cast<char*> (&gnfsMatrixRank), sizeof (long));

    if (!f.good()) {
        f.close();
        return false;
    }
    assert(columns == 0 || columns == gnfsMatrix.NumCols());

    freeVariables.SetLength(columns);

    if (columns == gnfsMatrix.NumCols()) {
        long wordsLength = (columns / NTL_BITS_PER_LONG) + 1;

        //Read in the wordsLength _ntl_ulongs before the last 2*sizeof(long) bytes. They are the freeVariables bit vector data...
        f.seekg(-(sizeof (_ntl_ulong) * wordsLength) - (2 * sizeof (long)), ios::end);

        for (long j = 0; j < wordsLength; ++j)
            f.read(reinterpret_cast<char*> (freeVariables.rep.elts() + j), sizeof (_ntl_ulong));

        if (!f.good()) {
            f.close();
            return false;
        }
    }

    f.close();

#ifdef GNFS_DEBUG
    assert(gnfsMatrixRank == gnfsMatrix.NumCols() - weight(freeVariables));
#endif

    return true;
}

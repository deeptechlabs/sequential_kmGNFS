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


#ifndef _GNFS_LINEARALGEBRA_H
#define	_GNFS_LINEARALGEBRA_H


#include <NTL/mat_GF2.h>

#include <GNFS/gnfs_global.h>
#include <GNFS/MatrixConstruction.h>
#include <GNFS/Relation.h>

class LinearAlgebra : public MatrixConstruction {
public:

    long gnfsMatrixRank;
    vec_GF2 freeVariables;

    static const char* reducedMatrixTextFile; //reduced matrix text file
    static const char* reducedMatrixBinaryFile; //reduced matrix binary file

    bool linearAlgebraInitialized;

public:

    LinearAlgebra()
    : MatrixConstruction(), linearAlgebraInitialized(false) {

    }

    bool init();
    bool reduceMatrixUsingGauss(bool readMatrixFromDisk = true, bool writeReducedMatrixToDisk = true, bool binaryIO = false);

protected:
    long gauss();
    bool getSolutions(mat_GF2& solutions, const long numberOfSolutions);
    bool getDependenciesFromSolutions(vec_vec_Relation& dependencies, const mat_GF2& solutions);
    bool getNextSolution(vec_GF2& solution);

    /* Gets a solution from the reduced matrix.
     * The values of the free variables should be set a priori in the solution vector
     */
    bool getSolution(vec_GF2& solution);

    bool writeMatrixText(const char* filename);
    bool writeMatrixBinary(const char* filename);
    bool readMatrixText(const char* filename);
    bool readMatrixBinary(const char* filename);



};

#endif	/* _GNFS_LINEARALGEBRA_H */


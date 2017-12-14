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


#include <iostream>
#include <string>
#include <cstdlib>
#include <climits>
#include <dirent.h>
#include <cerrno>
#include <cstring>

#include <NTL/mat_GF2.h>

#include <GNFS/Relation.h>
#include <GNFS/gnfs_utils.h>
#include <GNFS/gnfs_global.h>

#include <timer.h>

#include <GNFS/MatrixConstruction.h>

const char* MatrixConstruction::matrixTextFile = "IO/Matrix/constructed_matrix.txt";
const char* MatrixConstruction::matrixBinaryFile = "IO/Matrix/constructed_matrix.bin";

NTL_CLIENT


bool MatrixConstruction::init() {

    matrixConstructionInitialized = false;

    if (!Sieving::init())
        return false;


    if (!initializeRelationsDirectory())
        return false;

    totalNumberOfRelations = numOfFreeRelations + numOfSievingFoundRelations; //todo future
    cout << "Total number of relations: " << totalNumberOfRelations << endl;

    if (totalNumberOfRelations < numOfRequiredRelations) {
        cout << "The total number of found relations should be at least equal to (or larger than) the number of required relations (" << numOfRequiredRelations << ")..." <<
                endl << "Aborting..." << endl;
        return false;
    }

    kindOfReadRelations = 0; //start reading in full relations found during sieving

    matrixConstructionInitialized = true;
    return true;
}

bool MatrixConstruction::constructMatrix(bool binaryOutput) {

    if (!matrixConstructionInitialized) {
        cout << "Matrix construction step is not initialized properly! Please call init() first..." << endl;
        return false;
    }

    GF2 one;
    set(one);

    ZZ numerator;
    ZZ q;

    //if more of numOfRequiredRelations relations were found ( totalNumberOfRelations > numOfRequiredRelations )
    //then the excess relations are discarded
    long columns = numOfRequiredRelations;

    long rows = RFBSize + AFBSize + QCBSize + 1; //one row for the sign
    long signOffset = 0;
    long RFBOffset = 1;
    long AFBOffset = RFBSize + 1;
    long QCBOffset = RFBSize + AFBSize + 1;

    if (!monic) {
        ++rows; //one more row which will consist only of 1's and will ensure that the number of dependent relations is even
        ++signOffset;
        ++RFBOffset;
        ++AFBOffset;
        ++QCBOffset;
    }

    cout << "Number of rows: " << rows << endl;
    cout << "Number of columns: " << columns << endl;

    gnfsMatrix.kill();
    gnfsMatrix.SetDims(rows, columns);


    Relation rel;
    long length;

    for (long j = 0; j < columns; ++j) {
        if (!getNextRelationFromDisk(rel)) {
            return false;
        }

        //find the bit that corresponds to column j and create the appropriate mask for manipulating this bit
        long wj = j / NTL_BITS_PER_LONG;
        long bj = j - wj*NTL_BITS_PER_LONG;
        _ntl_ulong j_mask = 1UL << bj;

        if (!monic)
            gnfsMatrix[0].rep.elts()[wj] |= j_mask;


        if (rel.sign)
            gnfsMatrix[signOffset].rep.elts()[wj] |= j_mask;


        length = rel.rationalFactorsIndex.length();
        for (long i = 0; i < length; ++i) {
            if (rel.rationalFactors[i].b & 1)
                gnfsMatrix[rel.rationalFactorsIndex[i] + RFBOffset].rep.elts()[wj] |= j_mask;
        }

        length = rel.algebraicFactorsIndex.length();
        for (long i = 0; i < length; ++i) {
            if (rel.algebraicFactors[i].c & 1)
                gnfsMatrix[rel.algebraicFactorsIndex[i] + AFBOffset].rep.elts()[wj] |= j_mask;
        }


        if (rel.relPair.b == 0) {//handle free relations separately
            long p = rel.rationalFactors[0].a; //get corresponing p

            for (long i = 0; i < QCBSize; ++i) {

                conv(numerator, p);
                conv(q, QCB[i].a);
                rem(numerator, numerator, q); //reduce numerator modulo q

                if (Jacobi(numerator, q) == -1L)// check the Legendre symbol of (p / q)
                    gnfsMatrix[QCBOffset + i].rep.elts()[wj] |= j_mask;

            }
        }
        else {
            for (long i = 0; i < QCBSize; ++i) {

                conv(numerator, QCB[i].b); // s
                mul(numerator, numerator, rel.relPair.b); // bs
                sub(numerator, rel.relPair.a, numerator); // a - bs

                conv(q, QCB[i].a);
                rem(numerator, numerator, q); //reduce numerator modulo q

                if (Jacobi(numerator, q) == -1L)
                    gnfsMatrix[QCBOffset + i].rep.elts()[wj] |= j_mask;

            }
        }

        //cout << "Relation no " << i + 1 << " : " << rel << endl;

    }

    //cout << endl << gnfsMatrix << endl;

    finalizeMatrixConstruction();

    cout << "Writing constructed matrix to disk..." << endl;
    if (binaryOutput)
        return writeMatrixBinary(matrixBinaryFile);
    else
        return writeMatrixText(matrixTextFile);

    //return true;
}

bool MatrixConstruction::initializeRelationsDirectory() {

    if (pdir)
        rewinddir(pdir);
    else {
        pdir = opendir(relationsFileDirectoryPath.c_str());
        if (!pdir) {
            cout << "Error opening directory " << relationsFileDirectoryPath << endl;
            cout << "Make sure that the directory exists and is accessible..." << endl;
            return false;
        }
    }
    numOfReadRelations = 0;
    numOfReadRelationFiles = 0;
    kindOfReadRelations = 0;
    return true;
}

bool MatrixConstruction::getNextRelationFromDisk(Relation& rel) {

    if (!relationFileInStream.is_open()) {
        if (!openNextRelationsFile())
            return false;
    }


    relationFileInStream >> rel;
    relationFileInStream.ignore(INT_MAX, '\n');
    ++numOfReadRelations;

    SkipWhiteSpace(relationFileInStream);


    if (relationFileInStream.eof()) {
        relationFileInStream.close();
        relationFileInStream.clear();
    }

    if (!relationFileInStream.good()) {
        cout << "Error while reading the relations from disk..." << endl;
        return false;
    }

    return true;
}

bool MatrixConstruction::openNextRelationsFile() {

    switch (kindOfReadRelations) {

        case 0: //read in full relations found during sieving
            static size_t prefixLength = 0;

            if (prefixLength == 0) {
                prefixLength = strlen(relationsFilePrefix);
                cout << "Reading in full relations..." << endl;
            }

            errno = 0;
            while ((entry = readdir(pdir))) {
                if (strncmp(entry->d_name, relationsFilePrefix, prefixLength) == 0) {
                    currentRelationsFilePath = relationsFileDirectoryPath;
                    currentRelationsFilePath += entry->d_name;
                    relationFileInStream.open(currentRelationsFilePath.c_str());
                    if (!relationFileInStream.is_open()) {
                        cout << "Error opening " << currentRelationsFilePath << endl;
                        return false;
                    }
                    ++numOfReadRelationFiles;
                    return true;
                }
            }

            if (errno) {
                cout << "Error while reading the contents of " << relationsFileDirectoryPath << " ..." << endl;
                perror("The following error occured");
                return false;
            }

            //we have finished reading the full realations files

            //check to see whether all full relations have been read
            if (numOfReadRelationFiles != relationsFileCounter || numOfReadRelations != numOfSievingFoundRelations) {
                cout << "Could no read all full relations from disk. Some files and/or relations seem to be missing..." << endl;
                return false;
            }

            //switch to free relations file
            kindOfReadRelations = 1; //assigned only for the record, cause there is no actual use for this, since the execution continues to the next label anyway

        case 1: //read in free relations from free relations file

            cout << "Reading in free relations..." << endl;
            relationFileInStream.open(freeRelationsFile);
            if (!relationFileInStream.is_open()) {
                cout << "Error opening " << freeRelationsFile << endl;
                return false;
            }
            long tmpReadNumber;
            relationFileInStream >> tmpReadNumber; //read the number of free relations
            assert(numOfFreeRelations == tmpReadNumber);
            kindOfReadRelations = 2; //the next time read from the next kind of relations. Future use for partial relations...
            return true;

        case 2:
            return false;

    }

    return false;
}

void MatrixConstruction::finalizeMatrixConstruction() {

    closedir(pdir);
    pdir = 0;
    entry = 0;
}

bool MatrixConstruction::writeMatrixText(const char* filename) {

    ofstream f;
    f.open(filename);
    if (!f.is_open()) {
        cout << "Error opening " << filename << "..." << endl;
        return false;
    }

    f << gnfsMatrix << endl;

    if (f.good()) {
        f.close();
        return true;
    }

    cout << "Error while reading from " << filename << "..." << endl;
    f.close();
    return false;
}

bool MatrixConstruction::writeMatrixBinary(const char* filename) {

    long rows = gnfsMatrix.NumRows();
    long columns = gnfsMatrix.NumCols();
    if (rows == 0 || columns == 0) {
        cout << "Error while saving the matrix to disk. The matrix has zero rows and/or columns... Aborting..." << endl;
        return false;
    }

    long wordsLength = (columns / NTL_BITS_PER_LONG) + 1;

    ofstream f;
    f.open(filename, ios::out | ios::binary);
    if (!f.is_open()) {
        cout << "Error opening " << filename << "..." << endl;
        return false;
    }

    f.write(reinterpret_cast<char*> (&rows), sizeof (long));
    f.write(reinterpret_cast<char*> (&columns), sizeof (long));
    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < wordsLength; ++j)
            f.write(reinterpret_cast<char*> (gnfsMatrix[i].rep.elts() + j), sizeof (_ntl_ulong));
    }

    if (!f.good()) {
        f.close();
        return false;
    }

    f.close();
    return true;
}

bool MatrixConstruction::readMatrixText(const char* filename) {

    ifstream f;
    f.open(filename);
    if (!f.is_open()) {
        cout << "Error opening " << filename << "..." << endl;
        return false;
    }

    f >> gnfsMatrix;

    if (f.good()) {
        f.close();
        return true;
    }

    cout << "Error while reading from " << filename << "..." << endl;
    f.close();
    return false;
}

bool MatrixConstruction::readMatrixBinary(const char* filename) {

    ifstream f;
    f.open(filename, ios::in | ios::binary);
    if (!f.is_open()) {
        cout << "Error opening " << filename << "..." << endl;
        return false;
    }

    long rows = 0;
    long columns = 0;

    f.read(reinterpret_cast<char*> (&rows), sizeof (long));
    f.read(reinterpret_cast<char*> (&columns), sizeof (long));

    if (!f.good()) {
        f.close();
        return false;
    }

    gnfsMatrix.kill();
    gnfsMatrix.SetDims(rows, columns);
    long wordsLength = (columns / NTL_BITS_PER_LONG) + 1;

    for (long i = 0; i < rows; ++i) {
        for (long j = 0; j < wordsLength; ++j)
            f.read(reinterpret_cast<char*> (gnfsMatrix[i].rep.elts() + j), sizeof (_ntl_ulong));
    }

    if (!f.good()) {
        f.close();
        return false;
    }

    f.close();
    return true;
}

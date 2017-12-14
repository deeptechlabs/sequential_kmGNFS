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


#ifndef _GNFS_MATRIXCONSTRUCTION_H
#define	_GNFS_MATRIXCONSTRUCTION_H

#include <string>
#include <dirent.h>

#include <NTL/mat_GF2.h>

#include <GNFS/gnfs_global.h>
#include <GNFS/Sieving.h>
#include <GNFS/Relation.h>

class MatrixConstruction : public Sieving {
public:

    mat_GF2 gnfsMatrix;

    long totalNumberOfRelations;

    static const char* matrixTextFile; //constructed matrix text file
    static const char* matrixBinaryFile; //constructed matrix binary file


    bool matrixConstructionInitialized;

    DIR* pdir;
    struct dirent* entry;
    long kindOfReadRelations; //0 when reading in full relations found during sieving. 1 when reading in free relations

    long numOfReadRelations;
    long numOfReadRelationFiles;


public:

    MatrixConstruction()
    : Sieving(), matrixConstructionInitialized(false), pdir(0), entry(0) {

    }

    bool init();
    bool constructMatrix(bool binaryOutput = false);


protected:
    bool initializeRelationsDirectory();
    bool getNextRelationFromDisk(Relation& rel);
    bool openNextRelationsFile();
    void finalizeMatrixConstruction(); //close directory handles etc
    virtual bool writeMatrixText(const char* filename);
    virtual bool writeMatrixBinary(const char* filename);
    virtual bool readMatrixText(const char* filename);
    virtual bool readMatrixBinary(const char* filename);

};

#endif	/* _GNFS_MATRIXCONSTRUCTION_H */


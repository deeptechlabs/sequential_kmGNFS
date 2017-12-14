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
#include <cerrno>
#include <unistd.h>

#include <GNFS/gnfs_global.h>
#include <GNFS/FactorBases.h>

#include <timer.h>

using namespace std;


//The project name (GNFS_Factor_Bases_Setup) is a bit misleading because it constructs all the required bases,
//like QCB for example, and not only the factor bases

int main(int argc, char* argv[]) {

    if (argc != 2) {
        printGNFSUsage();
        return EXIT_FAILURE;
    }

    cout << "*** Setting current working dir to \"" << argv[1] << "\" ..." << endl;

    errno = 0;
    if (chdir(argv[1]) == -1) {
        cout << "Error while trying to set current working directory..." << endl;
        perror("The following error occured");
        printGNFSUsage();
        return EXIT_FAILURE;
    }

#ifdef GNFS_DEBUG
    cout << "$$$ Note: GNFS debug mode is on... $$$" << endl;
#endif

    cout << "*** Factor base construction step started..." << endl;
    cout << "*** Reading input params..." << endl;

    FactorBases fb;

    timer t;
    t.elapsed();

    cout << "*** Creating the Rational factor base (RFB)..." << endl;
    if (!fb.createRFB()) {
        cout << "*** An error occured while creating RFB... Exiting..." << endl;
        return EXIT_FAILURE;
    }


    cout << "*** Writing RFB to file(s)..." << endl;
    if (!fb.writeFBToFile(FactorBases::RFB_)) {
        cout << "*** An error occured while writing RFB to disk... Exiting..." << endl;
        return EXIT_FAILURE;
    }

    //    cout << "*** Reading RFB from disk..." << endl;
    //    if (!fb.readFBFromFile(FactorBases::RFB_)) {
    //        cout << "*** An error occured while reading RFB from disk... Exiting..." << endl;
    //        return EXIT_FAILURE;
    //    }



    cout << "*** Creating the Algebraic factor base (AFB)..." << endl;
    if (!fb.createAFB()) {
        cout << "*** An error occured while creating AFB... Exiting..." << endl;
        return EXIT_FAILURE;
    }

    cout << "*** Writing AFB to file(s)..." << endl;
    if (!fb.writeFBToFile(FactorBases::AFB_)) {
        cout << "*** An error occured while writing AFB to disk... Exiting..." << endl;
        return EXIT_FAILURE;
    }


    //    cout << "*** Reading AFB from disk..." << endl;
    //    if (!fb.readFBFromFile(FactorBases::AFB_)) {
    //        cout << "*** An error occured while reading AFB from disk... Exiting..." << endl;
    //        return EXIT_FAILURE;
    //    }


    cout << "*** Creating the Quadratic character base (QCB)..." << endl;
    if (!fb.createQCB()) {
        cout << "*** An error occured while creating QCB... Exiting..." << endl;
        return EXIT_FAILURE;
    }

    cout << "*** Writing QCB to file(s)..." << endl;
    if (!fb.writeFBToFile(FactorBases::QCB_)) {
        cout << "*** An error occured while writing QCB to disk... Exiting..." << endl;
        return EXIT_FAILURE;
    }


    //    cout << "*** Reading QCB from disk..." << endl;
    //    if (!fb.readFBFromFile(FactorBases::QCB_)) {
    //        cout << "*** An error occured while reading QCB from disk... Exiting..." << endl;
    //        return EXIT_FAILURE;
    //    }

    if (NTL_OVERFLOW(1, fb.RFBSize, fb.AFBSize)) {
        cout << "RFBSize + AFBSize shoud be smaller than or equal to NTL_OVFBND! Aborting..." << endl;
        return false;
    }

    long numOfRequiredRelations = fb.RFBSize + fb.AFBSize + 1; //one column for the sign
    if (!fb.monic)
        ++numOfRequiredRelations; //one column which will consist only of 1's and will ensure that the number of dependent relations is even

    if (NTL_OVERFLOW(1, numOfRequiredRelations, fb.QCBSize + fb.numOfExcessiveRelations)) {
        cout << "RFBSize + AFBSize + QCBSize + numOfExcessRelations shoud be smaller than or equal to NTL_OVFBND! Aborting..." << endl;
        return false;
    }

    numOfRequiredRelations += fb.QCBSize + fb.numOfExcessiveRelations;

    cout << "\nRequired relations to be found: " << numOfRequiredRelations << endl;

    double ti = t.elapsed();
    cout << "Total Time: " << ti << endl;

    cout << endl << "*** All done! ***" << endl;

#ifdef GNFS_DEBUG
    cout << "$$$ Note: GNFS debug mode was on... $$$" << endl;
#endif


    //system("PAUSE");
    return EXIT_SUCCESS;
}

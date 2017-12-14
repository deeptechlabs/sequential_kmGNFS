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
#include <GNFS/Initialization.h>

#include <timer.h>

using namespace std;

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

    cout << "*** Initialization step started..." << endl;
    cout << "*** Reading input params..." << endl;

    Initialization init;

    cout << "*** Running some initial tests..." << endl;

    timer t;
    t.elapsed();

    if (init.runInitialTests())
        cout << "*** All tests passed..." << endl;
    else {
        cout << "*** Some tests failed... Aborting..." << endl;
        return EXIT_FAILURE;

    }

    if (init.performPolynomialSelection) {
        cout << "*** Performing polynomial selection..." << endl;

        if (!init.polynomialSelectionStep()) {
            cout << "An error occurred while trying to find the best rated polynomial through an extensive search..." << endl;
            return EXIT_FAILURE;
        }

        cout << "*** Polynomial selection completed..." << endl;
    }
    else
        cout << "*** Skipping polynomial selection step..." << endl;



    cout << "*** Writing output parameters file..." << endl;
    init.writeParamsFile();
    cout << endl << "*** All done! ***" << endl;

#ifdef GNFS_DEBUG
    cout << "$$$ Note: GNFS debug mode was on... $$$" << endl;
#endif


    double time = t.elapsed();
    cout << endl << "Total Time: " << time << endl;

    //system("PAUSE");
    return EXIT_SUCCESS;
}

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
#include <GNFS/LineSieving.h>

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

    cout << "*** Sieving step started..." << endl;
    cout << "*** Reading input params..." << endl;

    LineSieving ls;

    timer t;
    t.elapsed();

    cout << "*** Initializing sieving step..." << endl;
    if (!ls.init()) {
        cout << "*** An error occured while initializing sieving step... Exiting..." << endl;
        return EXIT_FAILURE;
    }



    cout << "*** Executing the sieving step..." << endl;
    if (!ls.doLineSieving()) {
        cout << "*** An error occured while executing the sieving step... Exiting..." << endl;
        return EXIT_FAILURE;
    }

    double ti = t.elapsed();
    cout << "Total Time: " << ti << endl;


    cout << endl << "*** All done! ***" << endl;

#ifdef GNFS_DEBUG
    cout << "$$$ Note: GNFS debug mode was on... $$$" << endl;
#endif


    //system("PAUSE");
    return EXIT_SUCCESS;
}

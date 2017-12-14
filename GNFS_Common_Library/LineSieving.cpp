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


#include <GNFS/LineSieving.h>

#include <timer.h>

const char* LineSieving::relationsFileInfix = "_start_b_";

static int terminateSieving = 0;

NTL_CLIENT

/* This signal handler will be used to catch both SIGINT and SIGUSR1 signals.
 * Its purpose is to notify the sieving algorithm that the user has requested
 * termination of the sieving process, and let the sieving algorithm finish
 * sieving the current b value, save its status and then exit normally.
 */
static void catch_signal(int signum) {

    terminateSieving++;
    if (terminateSieving > 2) {
        cout << "Aborting immediately..." << endl;
        exit(1);
    }
    cout << "Terminating gracefully... Please wait..." << endl;
    signal(SIGINT, catch_signal);
    signal(SIGUSR1, catch_signal);
}

bool LineSieving::init() {

    bool discardPreviousSievingResults = false;

#ifdef GNFS_DEBUG
    cout << "### Warning: Debug mode is ON. Previous sieving results will be erased... ###" << endl;
    discardPreviousSievingResults = true;
#endif

    lineSevingInitialized = false;
    if (!Sieving::init(discardPreviousSievingResults))
        return false;

    cout << "Last b value that has been previously sieved is " << last_b << "..." << endl;

    lineSevingInitialized = true;
    return true;
}

bool LineSieving::doLineSievingWithLogarithms() {

    if (!lineSevingInitialized) {
        cout << "Sieving step is not initialized properly! Please call init() first..." << endl;
        return false;
    }

    signal(SIGINT, catch_signal);
    signal(SIGUSR1, catch_signal);

    /*timer t;
    t.elapsed();*/

    //set up rational sieve//

    sieve_R.FB = &RFB;
    sieve_R.FBSize = RFBSize;
    sieve_R.whichBase = RFB_;

    sieve_R.b_x_r_mod_p = new (nothrow) int_type[RFBSize];
    if (sieve_R.b_x_r_mod_p == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    int_type* b_x_m_mod_p = sieve_R.b_x_r_mod_p;


    sieve_R.a_min_mod_p = new (nothrow) int_type[RFBSize];
    if (sieve_R.a_min_mod_p == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    int_type* a_min_mod_p_R = sieve_R.a_min_mod_p;


    sieve_R.logP = new (nothrow) sieve_type[RFBSize];
    if (sieve_R.logP == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    sieve_type* logP_R = sieve_R.logP;


    long a_min = -max_a;
    long rowLength = 2 * max_a + 1;
    long b = max(min_b, last_b + 1);


    for (long i = 0; i < RFBSize; ++i) {
        long p = RFB[i].a;
        b_x_m_mod_p[i] = MulMod(b % p, RFB[i].b, p); // RFB[i].b == m % p

        a_min_mod_p_R[i] = reduce(a_min, p);
    }

    sieve_R.row = new (nothrow) sieve_type[rowLength];
    if (sieve_R.row == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    sieve_type *row_R = sieve_R.row;
    sieve_R.rowLength = rowLength;

    sieve_R.a_min = -max_a;
    sieve_R.a_max = max_a;
    sieve_R.b_min = b;
    sieve_R.b_max = max_b;

    ZZ b_x_m; //b*m
    mul(b_x_m, m, b);
    ZZ F2_1, F2_2; // a +/- b*m


    //set up algebraic sieve//

    sieve_A.FB = &AFB;
    sieve_A.FBSize = AFBSize;
    sieve_A.whichBase = AFB_;

    sieve_A.b_x_r_mod_p = new (nothrow) int_type[AFBSize];
    if (sieve_A.b_x_r_mod_p == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    int_type* b_x_r_mod_p = sieve_A.b_x_r_mod_p;

    sieve_A.a_min_mod_p = new (nothrow) int_type[AFBSize];
    if (sieve_A.a_min_mod_p == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    int_type* a_min_mod_p_A = sieve_A.a_min_mod_p;

    sieve_A.logP = new (nothrow) sieve_type[AFBSize];
    if (sieve_A.logP == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    sieve_type* logP_A = sieve_A.logP;

    for (long i = 0; i < AFBSize; ++i) {
        long p = AFB[i].a;
        if (p == AFB[i].b) {
            //handle projective root entries     
            b_x_r_mod_p[i] = 0;
        }
        else
            b_x_r_mod_p[i] = MulMod(b % p, AFB[i].b, p); // AFB[i].b == r % p

        a_min_mod_p_A[i] = reduce(a_min, p);
    }

    sieve_A.row = new (nothrow) sieve_type[rowLength];
    if (sieve_A.row == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    sieve_type *row_A = sieve_A.row;
    sieve_A.rowLength = rowLength;

    sieve_A.a_min = -max_a;
    sieve_A.a_max = max_a;
    sieve_A.b_min = b;
    sieve_A.b_max = max_b;

    ZZ F1_1, F1_2; //homogeneous F

    //compute log bases and logs of primes for both sieves, if necessary//
    long tmp;
    if ((tmp = b % logBaseUpdateRate) != 0) {
        //if b % logBaseUpdateRate == 0 the update will take place right below in the main loop
        computeLogarithms(sieve_R, b + logBaseUpdateRate - tmp - 1);
        computeLogarithms(sieve_A, b + logBaseUpdateRate - tmp - 1);
    }



    //skip small primes for the rational sieve//
    long primeOffset_R = 0;
    while (RFB[primeOffset_R].a < ignoreSmallPrimesBound)
        ++primeOffset_R;
    sieve_R.startingPrimeOffset = primeOffset_R;

    //skip small primes for the algebraic sieve//
    long primeOffset_A = 0;
    while (AFB[primeOffset_A].a < ignoreSmallPrimesBound)
        ++primeOffset_A;
    sieve_A.startingPrimeOffset = primeOffset_A;


    //sieve by rows//

    Relation tmpRel;
    vec_pair_intType_intType reports;


    for (/*b = max(min_b, last_b + 1)*/; b <= max_b; ++b) {

        if (terminateSieving)
            break;

        //recompute log bases and logs of primes for both sieves, if necessary//
        if (b % logBaseUpdateRate == 0) {
            computeLogarithms(sieve_R, b + logBaseUpdateRate - 1);
            computeLogarithms(sieve_A, b + logBaseUpdateRate - 1);
        }


        //prepare rational row//

        bzero(row_R, rowLength * sizeof (sieve_type));


        //sieve the rational row//

        long p;
        sieve_type logP;
        long startingOffset;

        for (long i = primeOffset_R; i < RFBSize; ++i) {
            p = RFB[i].a;
            logP = logP_R[i];

            //cout << p << endl;
            //long tmp = a_min % p;
            /* if (tmp < 0)
                 tmp += p;*/
            startingOffset = calculateSievingOffset(a_min_mod_p_R[i], b_x_m_mod_p[i], p);


#ifdef GNFS_DEBUG
            assert(startingOffset >= 0);
            assert(startingOffset < p);
#endif




            for (long off = startingOffset; off < rowLength; off += p) {
                row_R[off] += logP;
            }

        }




        //prepare algebraic row//

        /* find the projective roots that divide b, and compute
           their combined log value. This is added to the cutoff
           for the entire sieve later */

        long commonFactors = f.rep[d] % b; // f.rep[d] is a bit faster than coeff(f, d)
        commonFactors = GCD(b, commonFactors);
        sieve_A.projectiveRootsLogBias = static_cast<sieve_type> (log(commonFactors, sieve_A.logOfBase));


        memset(row_A, sieve_A.projectiveRootsLogBias, rowLength * sizeof (sieve_type));
        //bzero(row_A, rowLength * sizeof (sieve_type));

        //sieve the algebraic row//

        for (long i = primeOffset_A; i < AFBSize; ++i) {
            p = AFB[i].a;
            logP = logP_A[i];


            //skip if p divides b or if the entry is a projective root//
            if (b_x_r_mod_p[i] == 0 && (AFB[i].b == p || b % p == 0))
                continue;



            startingOffset = calculateSievingOffset(a_min_mod_p_A[i], b_x_r_mod_p[i], p);

            //cout << "prime: " << p << " log: " << (unsigned int)logP << " st_off: " << startingOffset << endl;

#ifdef GNFS_DEBUG
            assert(startingOffset >= 0);
            assert(startingOffset < p);
            //assert(divide(row_A[startingOffset], p));
            /*if(!divide(row_A[startingOffset], p))
            {
                cout << "oooppes" << endl;
                cout << startingOffset << endl;
                cout << p << endl;
                cout << AFB[i] << endl;
                cout << b << endl;
                cout << evalHomogeneousF(f,startingOffset + a_min,b) << endl;
                cout << startingOffset + a_min << endl;
                cout << a_min_mod_p_A[i] << endl;
                cout << b_x_r_mod_p[i] << endl;
                cout << row_A[startingOffset] << endl;
                cout << (startingOffset + a_min) % p << endl;
                exit(0);
            }*/
#endif


            for (long off = startingOffset; off < rowLength; off += p) {
                row_A[off] += logP;
            }



        }

        /*cout << "---------RFB logs-----------" << endl;
        for (long i = 0; i < RFBSize; ++i) {
            cout << "prime: " << RFB[i].a << " log: "  << (int)logP_R[i] << endl;
        }
        cout << "--------AFB logs----------" << endl;
        for (long i = 0; i < AFBSize; ++i) {
            cout << "prime: " << AFB[i].a << " log: "  << (int)logP_A[i] << endl;
        }*/
        /*for (long i = 0, a = a_min; i < rowLength; ++i, ++a) {
            cout << (int)row_R[i] << endl;
        }*/
        /*for (long i = 0, a = a_min; i < rowLength; ++i, ++a) {
            cout << i << ":  " << (int)row_A[i] << endl;
        }*/

        if (terminateSieving)
            break;

        //scan for reports//

        long cutoff_R;
        long cutoff_A;
        long a = a_min;
        long length, innerLength;
        long right, innerRight;

        for (long left = 0; left < rowLength; left += cuttofUpdateRate_R) {
            right = left + cuttofUpdateRate_R - 1;
            if (right >= rowLength)
                right = rowLength - 1;
            length = right - left + 1;


            //sub(F2_1, a, b_x_m);
            //sub(F2_2, a + length - 1, b_x_m);
            F2_1 = F2_2 = b_x_m;

            abs(F2_1, F2_1);
            abs(F2_2, F2_2);

            //this should occur only when factoring really small numbers...
            if (IsZero(F2_1)) ++F2_1;
            if (IsZero(F2_2)) ++F2_2;


            cutoff_R = (log(F2_1, sieve_R.logOfBase) + log(F2_2, sieve_R.logOfBase)) / 2;
            cutoff_R -= sieve_R.fudgeFactor;

            for (long innerLeft = left; innerLeft <= right; innerLeft += cuttofUpdateRate_A) {
                innerRight = innerLeft + cuttofUpdateRate_A - 1;
                if (innerRight > right)
                    innerRight = right;
                innerLength = innerRight - innerLeft + 1;

                evalHomogeneousF(F1_1, f, a, b);
                evalHomogeneousF(F1_2, f, a + innerLength - 1, b);

                abs(F1_1, F1_1);
                abs(F1_2, F1_2);

                cutoff_A = (log(F1_1, sieve_A.logOfBase) + log(F1_2, sieve_A.logOfBase)) / 2;
                cutoff_A -= sieve_A.fudgeFactor;

                for (long i = innerLeft; i <= innerRight; ++i, ++a) {
                    if (row_A[i] <= cutoff_A || row_R[i] <= cutoff_R || GCD(a, b) != 1)
                        continue;

                    ++numOfReports;
                    tmpRel.relPair.a = a;
                    tmpRel.relPair.b = b;
                    append(reports, tmpRel.relPair);

                    /*if (trialDivideRelation(tmpRel, b_x_m, b_x_m_mod_p, b_x_r_mod_p)) {
                        if (!writeRelationToFile(tmpRel))
                            return false;
                         ++numOfSievingFoundRelations;
                    }*/
                }
            }

        }

        //trial divide reports to find smooth pairs//
        long reportsLength = reports.length();
        for (long i = 0; i < reportsLength; ++i) {
            tmpRel.relPair = reports[i];

            if (trialDivideRelation(tmpRel, b_x_m, b_x_m_mod_p, b_x_r_mod_p)) {
                if (!writeRelationToFile(tmpRel))
                    return false;
                ++numOfSievingFoundRelations;
            }
        }
        reports.QuickSetLength(0);


        //this is for the the next algebraic row//
        for (long i = 0; i < AFBSize; ++i) {
            b_x_r_mod_p[i] = AddMod(b_x_r_mod_p[i], AFB[i].b, AFB[i].a);

        }

        //this is for the the next rational row//
        for (long i = 0; i < RFBSize; ++i) {
            b_x_m_mod_p[i] = AddMod(b_x_m_mod_p[i], RFB[i].b, RFB[i].a);

        }
        add(b_x_m, b_x_m, m);

        if (b % 50 == 0) {
            cout << "Sieved for b " << b << " of " << max_b << ". Found " << numOfSievingFoundRelations << " full relations (" << numOfReports << " reports)..." << endl;
        }

    }

    last_b = b - 1;

    /*double ti = t.elapsed();
    cout << "Time: " << ti << endl;*/

    cout << "Total number of reports: " << numOfReports << endl;
    cout << "Found a total of " << numOfSievingFoundRelations << " full relations through line sieving with logarithms..." << endl;
    cout << "Total gathered relations: " << numOfFreeRelations + numOfSievingFoundRelations << endl;


    return finalizeSieving();


}

bool LineSieving::doSimpleLineSieving() {

    if (!lineSevingInitialized) {
        cout << "Sieving step is not initialized properly! Please call init() first..." << endl;
        return false;
    }

    signal(SIGINT, catch_signal);
    signal(SIGUSR1, catch_signal);
    /*timer t;
    t.elapsed();*/

    //set up rational sieve//

    sieve_R.b_x_r_mod_p = new (nothrow) int_type[RFBSize];
    if (sieve_R.b_x_r_mod_p == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    int_type* b_x_m_mod_p = sieve_R.b_x_r_mod_p;


    sieve_R.a_min_mod_p = new (nothrow) int_type[RFBSize];
    if (sieve_R.a_min_mod_p == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    int_type* a_min_mod_p_R = sieve_R.a_min_mod_p;


    long a_min = -max_a;
    long rowLength = 2 * max_a + 1;
    long b = max(min_b, last_b + 1);

    for (long i = 0; i < RFBSize; ++i) {
        long p = RFB[i].a;
        b_x_m_mod_p[i] = MulMod(b % p, RFB[i].b, p); // RFB[i].b == m % p

        a_min_mod_p_R[i] = reduce(a_min, p);
    }


    ZZ b_x_m; //b*m
    mul(b_x_m, m, b);
    ZZ F2; // a +/- b*m

    ZZ largestCell;
    mul(largestCell, m, max_b);
    add(largestCell, largestCell, max_a);

    ZZVec row_R;
    row_R.SetSize(rowLength, largestCell.size());


    //set up algebraic sieve//

    sieve_A.b_x_r_mod_p = new (nothrow) int_type[AFBSize];
    if (sieve_A.b_x_r_mod_p == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    int_type* b_x_r_mod_p = sieve_A.b_x_r_mod_p;

    sieve_A.a_min_mod_p = new (nothrow) int_type[AFBSize];
    if (sieve_A.a_min_mod_p == 0) {
        cout << "Memory allocation error!!" << endl;
        return false;
    }
    int_type* a_min_mod_p_A = sieve_A.a_min_mod_p;

    for (long i = 0; i < AFBSize; ++i) {
        long p = AFB[i].a;
        if (p == AFB[i].b) {
            //handle projective root entries
            b_x_r_mod_p[i] = 0;
        }
        else
            b_x_r_mod_p[i] = MulMod(b % p, AFB[i].b, p); // AFB[i].b == r % p

        a_min_mod_p_A[i] = reduce(a_min, p);
    }

    vec_ZZ row_A;
    row_A.SetLength(rowLength);


    Relation tmpRel;

    //sieve by rows//

    for (/*b = max(min_b, last_b + 1)*/; b <= max_b; ++b) {

        if (terminateSieving)
            break;

        //prepare rational row//

        sub(F2, a_min, b_x_m);

        for (long i = 0; i < rowLength; ++i) {
            row_R[i] = F2;
            ++F2;
        }


        //sieve the rational row//

        long p;
        long startingOffset;
        for (long i = 0; i < RFBSize; ++i) {
            p = RFB[i].a;
            //cout << p << endl;
            //long tmp = a_min % p;
            /* if (tmp < 0)
                 tmp += p;*/
            startingOffset = calculateSievingOffset(a_min_mod_p_R[i], b_x_m_mod_p[i], p);


#ifdef GNFS_DEBUG
            assert(startingOffset >= 0);
            assert(startingOffset < p);
            assert(divide(row_R[startingOffset], p));
#endif




            for (long i = startingOffset; i < rowLength; i += p) {
                if (IsZero(row_R[i]))
                    continue;
                do {
                    div(row_R[i], row_R[i], p);
                } while (row_R[i] % p == 0);
            }

        }

        //prepare algebraic row//

        for (long i = 0, a = a_min; i < rowLength; ++i, ++a) {
            evalHomogeneousF(row_A[i], f, a, b);
        }


        //sieve the algebraic row//

        for (long i = 0; i < AFBSize; ++i) {
            p = AFB[i].a;



            //skip if p divides b or if the entry is a projective root// 
            if (b_x_r_mod_p[i] == 0 && (AFB[i].b == p || b % p == 0))
                continue;

            //CAUTION: we skip all projective roots, so we miss ALL entries that are divided by them
            //(this happens for the entire row when p divides b)

            startingOffset = calculateSievingOffset(a_min_mod_p_A[i], b_x_r_mod_p[i], p);


#ifdef GNFS_DEBUG
            assert(startingOffset >= 0);
            assert(startingOffset < p);
            assert(divide(row_A[startingOffset], p));
            /* if (!divide(row_A[startingOffset], p)) {
                 cout << "oooppes" << endl;
                 cout << startingOffset << endl;
                 cout << p << endl;
                 cout << AFB[i] << endl;
                 cout << b << endl;
                 ZZ tmp;
                 evalHomogeneousF(tmp, f, startingOffset + a_min, b);
                 cout << tmp << endl;
                 cout << startingOffset + a_min << endl;
                 cout << a_min_mod_p_A[i] << endl;
                 cout << b_x_r_mod_p[i] << endl;
                 cout << row_A[startingOffset] << endl;
                 cout << row_A[startingOffset + 1] << endl;
                 cout << row_A[startingOffset + 2] << endl;
                 cout << (startingOffset + a_min) % p << endl;
                 exit(0);
             }*/
#endif


            for (long i = startingOffset; i < rowLength; i += p) {
                if (IsZero(row_A[i]))
                    continue;
#ifdef GNFS_DEBUG
                assert(divide(row_A[i], p));
#endif
                while (row_A[i] % p == 0)
                    div(row_A[i], row_A[i], p);
            }



        }


        if (terminateSieving)
            break;

        //check for smooth pairs//

        for (long i = 0, a = a_min; i < rowLength; ++i, ++a) {

            if (IsOne(abs(row_R[i])) && IsOne(abs(row_A[i])) && GCD(a, b) == 1) {
                tmpRel.relPair.a = a;
                tmpRel.relPair.b = b;

#ifdef GNFS_DEBUG
                assert(trialDivideRelation(tmpRel, b_x_m, b_x_m_mod_p, b_x_r_mod_p));
#else
                trialDivideRelation(tmpRel, b_x_m, b_x_m_mod_p, b_x_r_mod_p);
#endif


                //cout << "Relation: " << a << "  " << b << endl;
                if (!writeRelationToFile(tmpRel))
                    return false;
                ++numOfSievingFoundRelations;
            }
        }


        //this is for the the next algebraic row//
        for (long i = 0; i < AFBSize; ++i) {
            b_x_r_mod_p[i] = AddMod(b_x_r_mod_p[i], AFB[i].b, AFB[i].a);

        }

        //this is for the the next rational row//
        for (long i = 0; i < RFBSize; ++i) {
            b_x_m_mod_p[i] = AddMod(b_x_m_mod_p[i], RFB[i].b, RFB[i].a);

        }
        add(b_x_m, b_x_m, m);

        if (b % 50 == 0) {
            cout << "Sieved for b: " << b << " of total " << max_b << ". Found " << numOfSievingFoundRelations << " full relations until now..." << endl;
        }


    }

    last_b = b - 1;

    /*double ti = t.elapsed();
    cout << "Time: " << ti << endl;*/

    cout << "Found a total of " << numOfSievingFoundRelations << "  full relations through line sieving with trial division..." << endl;


    return finalizeSieving();
}

bool LineSieving::trialDivideRelation(Relation& rel, const ZZ& b_x_m, const int_type* b_x_r_mod_p_R, const int_type* b_x_r_mod_p_A) { //Caution: a,b should be coprime!

    static ZZ F1, F2;
    static pair_intType_intType pair;
    static triple_intType triple;

    long a = rel.relPair.a;
    long b = rel.relPair.b;

    rel.softReset();


    sub(F2, a, b_x_m);


    //this should occur only when factoring really small numbers...
    if (IsZero(F2))
        return false;

    switch (sign(F2)) {
        case -1 :
            rel.sign = true;
            break;
        default:
            rel.sign = false;
    }

    abs(F2, F2);

    int_type a_mop_p;
    long p = 0;
    long exponentCounter = 0;
    for (long i = 0; i < RFBSize; ++i) {

        //if (p != RFB[i].a) { //useless yet
        p = RFB[i].a;
        a_mop_p = reduce(a, p);
        //}


        /*if (p == RFB[i].b) { //useless yet
            if (b % p != 0) 
                continue;
            
        }
        else*/ if (a_mop_p != b_x_r_mod_p_R[i])
            continue;


        while (F2 % p == 0) {
            div(F2, F2, p);
            ++exponentCounter;
        }

#ifdef GNFS_DEBUG
        assert(exponentCounter > 0);
#endif


        pair.a = p;
        pair.b = exponentCounter;
        rel.appendRationalFactor(pair, i);
        exponentCounter = 0;

        if (IsOne(F2))
            break;
    }

    if (!IsOne(F2))
        return false;


    evalHomogeneousF(F1, f, a, b);
    abs(F1, F1);


    p = 0;
    for (long i = 0; i < AFBSize;) {

        if (p != AFB[i].a) {
            p = AFB[i].a;
            a_mop_p = reduce(a, p);
        }

        /*if (b_x_r_mod_p_A[i] == 0 && (AFB[i].b == p || b % p == 0)) {
            ++i;
            continue;
        }*/

        /*  if (((p != AFB[i].b && a_mop_p == b_x_r_mod_p_A[i]) || (p == AFB[i].b && b % p == 0)) && F1 % p != 0) {
              cout << "fuck!!!!!" << endl;
              exit(1);
          }*/

        if (p == AFB[i].b) {

            if (b % p != 0) {
                ++i;
                continue;
            }
        }
        else if (a_mop_p != b_x_r_mod_p_A[i]) {
            ++i;
            continue;
        }

        /*if ((a_mop_p != b_x_r_mod_p_A[i]) && (p == AFB[i].b && b % p != 0)) {
            ++i;
            if(i < AFBSize && AFB[i].a != p && F1 % p == 0) {
                cout << p << endl;
                if(p == AFB[i-1].b && b % p != 0)
                    cout << "one" << endl;
                else
                    cout << "two" << endl;

                cout << AFB[i-1] << endl;
                cout << rel.relPair << endl;
                cout << a_mop_p << endl;

                assert(1==2);
            }
            continue;

        }*/


        while (F1 % p == 0) {

            div(F1, F1, p);
            ++exponentCounter;

        }

#ifdef GNFS_DEBUG
        assert(exponentCounter > 0);
#endif

        // if (exponentCounter > 0) {

        //if p is a projective root skip the step of finding which (p,r) is responsible for dividing F1
        /*           if (p != AFB[i].b || b % p != 0) {
                       if (p == AFB[i].b) {
                           cout << "fkljflsd" << endl;
                           ++i;
       #ifdef GNFS_DEBUG
                           assert(p == AFB[i].a);
       #endif
                       }*/

        /*if (b % p == 0)
            assert(1 == 2);*/

        //a_mop_p = reduce(a, p);
        /*
     
                     while (a_mop_p != b_x_r_mod_p_A[i])
    
                     {
                         ++i;
     #ifdef GNFS_DEBUG
                         assert(i < AFBSize);
                         assert(p == AFB[i].a);
                   
     #endif
                     }*/
        //      }

        /* else if (p == AFB[i].b)
             assert(b % p == 0);*/


        triple.a = p;
        triple.b = AFB[i].b;
        triple.c = exponentCounter;
        rel.appendAlgebraicFactor(triple, i);
        exponentCounter = 0;

        if (IsOne(F1))
            break;

        //    }

        do {
            ++i;
            /* if (i == AFBSize)
                 continue;
             long pp = AFB[i].a;
             
             if (p == pp && a_mop_p == b_x_r_mod_p_A[i]) {
                 cout << rel.relPair << endl;
                 cout << AFB[i].b << endl;
                 cout << a_mop_p << "  " << b_x_r_mod_p_A[i] << endl;
                 cout << p << "  " << flag << endl;
                 ZZ norm;
                 evalHomogeneousF(norm, f, a, b);
                 cout << norm << endl;
                 if (norm % p == 0)
                     cout << "divides!!" << endl;
                 assert(1 == 2);
             }*/
        } while (p == AFB[i].a && i < AFBSize);

    }

    if (!IsOne(F1))
        return false;

    return true;
}

bool LineSieving::writeRelationToFile(const Relation& rel) {

    static long relationsPerFileCounter = 0;

    if (!relationFileOutStream.is_open()) {

        ostringstream osstr;
        osstr << relationsFileDirectoryPath << relationsFilePrefix << ++relationsFileCounter << relationsFileInfix << rel.relPair.b << relationsFileSuffix;

        currentRelationsFilePath = osstr.str();

        relationFileOutStream.clear();
        relationFileOutStream.open(currentRelationsFilePath.c_str());
        if (!relationFileOutStream.is_open()) {
            cout << "Error opening " << currentRelationsFilePath << endl;
            cout << "Make sure that the directory " << relationsFileDirectoryPath << " exists and is accessible..." << endl;
            return false;
        }

        relationsPerFileCounter = 0;

    }
    else if (relationsPerFileCounter == relationsPerFileBound) {

        relationFileOutStream.close();
        if (!updateSievingStatus())
            return false;

        ostringstream osstr;
        osstr << relationsFileDirectoryPath << relationsFilePrefix << ++relationsFileCounter << relationsFileInfix << rel.relPair.b << relationsFileSuffix;

        currentRelationsFilePath = osstr.str();

        relationFileOutStream.clear();
        relationFileOutStream.open(currentRelationsFilePath.c_str());
        if (!relationFileOutStream.is_open()) {
            cout << "Error opening " << currentRelationsFilePath << endl;
            cout << "Make sure that the directory " << relationsFileDirectoryPath << " exists and is accessible..." << endl;
            return false;
        }


        relationsPerFileCounter = 0;
    }

    relationFileOutStream << rel << '\n';

    ++relationsPerFileCounter;

    if (!relationFileOutStream.good()) {
        cout << "Error while writing to " << currentRelationsFilePath << " ..." << endl;
        return false;
    }


    return true;
}

bool LineSieving::getSievingStatus() {

    numOfSievingFoundRelations = 0;
    relationsFileCounter = 0;
    numOfReports = 0;
    last_b = 0;

    sievingStatusFileInStream.open(sievingStatusFilePath.c_str());
    if (!sievingStatusFileInStream.is_open()) {//we assume that the file does not exist, thus there are not any previous results
        sievingStatusFileInStream.clear();
        return true;
    }

    sievingStatusFileInStream.peek();
    if (sievingStatusFileInStream.eof()) {
        sievingStatusFileInStream.close();
        return true;
    }

    sievingStatusFileInStream >> numOfSievingFoundRelations;
    sievingStatusFileInStream.ignore(INT_MAX, '\n');
    sievingStatusFileInStream >> relationsFileCounter;
    sievingStatusFileInStream.ignore(INT_MAX, '\n');
    sievingStatusFileInStream >> numOfReports;
    sievingStatusFileInStream.ignore(INT_MAX, '\n');
    sievingStatusFileInStream >> last_b;
    sievingStatusFileInStream.ignore(INT_MAX, '\n');

    if (!sievingStatusFileInStream.good()) {
        cout << "Error while reading from " << sievingStatusFilePath << " ..." << endl;
        sievingStatusFileInStream.close();
        return false;
    }

    sievingStatusFileInStream.close();

    return true;
}

bool LineSieving::updateSievingStatus() {

    if (!Sieving::updateSievingStatus())
        return false;

    if (!sievingStatusFileOutStream.is_open()) {

        sievingStatusFileOutStream.clear();
        sievingStatusFileOutStream.open(sievingStatusFilePath.c_str(), ios::app);
        if (!sievingStatusFileOutStream.is_open()) {
            cout << "Error opening " << sievingStatusFilePath << endl;
            cout << "Make sure that the directory " << relationsFileDirectoryPath << " exists and is accessible..." << endl;
            return false;
        }
    }

    sievingStatusFileOutStream.seekp(0, ios::end);

    sievingStatusFileOutStream << last_b << "\t# Last b that has been previously sieved" << endl;

    if (!sievingStatusFileOutStream.good()) {
        cout << "Error while writing to " << sievingStatusFilePath << " ..." << endl;
        return false;
    }


    return true;
}

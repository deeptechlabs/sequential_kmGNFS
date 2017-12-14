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


#include <cmath>
#include <vector>
#include <queue>
#include <iostream>
#include <fstream>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include <GNFS/gnfs_utils.h>
#include <GNFS/PolynomialSelection.h>

NTL_CLIENT

        const char* PolynomialSelection::bestPolynomialsFile = "IO/Initialization/best_polynomials_found.txt"; //The file which contains the best polynomials found through the polynomial selection process

ZZX PolynomialSelection::findMonicPolynomial(const ZZ& m) const {

    ZZX poly(INIT_SIZE, d + 1);

    ZZ power_of_m(INIT_SIZE, n.size());
    ZZ c_i(INIT_SIZE, n.size());
    ZZ tmp_n = n;

    power(power_of_m, m, d);
    sub(tmp_n, tmp_n, power_of_m);
    SetCoeff(poly, d, 1);

    for (long i = d - 1; i >= 0; --i) {
        div(power_of_m, power_of_m, m); //m^i = m^(i+1) / m
        DivRem(c_i, tmp_n, tmp_n, power_of_m); //c_i = floor(tmp_n / power_of_m); tmp_n = tmp_n - c_i*power_of_m;
        SetCoeff(poly, i, c_i);
    }

    return poly;
}

//todo

ZZX PolynomialSelection::findNonSkewedNonMonicPolynomial(const ZZ& m) const {

    //todo
    //
    //Thelei ftiaksimo

    ZZX poly(INIT_SIZE, d + 1);

    ZZ power_of_m(INIT_SIZE, n.size());
    ZZ c_i(INIT_SIZE, n.size());
    ZZ tmp_n = n;

    power(power_of_m, m, d);
    DivRem(c_i, tmp_n, tmp_n, power_of_m); //c_i = floor(tmp_n / power_of_m); tmp_n = tmp_n - c_i*power_of_m;
    SetCoeff(poly, d, c_i);

    for (long i = d - 1; i >= 0; --i) {
        div(power_of_m, power_of_m, m); //m^i = m^(i+1) / m
        DivRem(c_i, tmp_n, tmp_n, power_of_m); //c_i = floor(tmp_n / power_of_m); tmp_n = tmp_n - c_i*power_of_m;
        SetCoeff(poly, i, c_i);
    }

    return poly;
}

ZZX PolynomialSelection::findNonSkewedNonMonicPolynomialWithLeadingCoefKnown(const ZZ& m, const ZZ &leadingCoef) const {

    ZZX poly(INIT_SIZE, d + 1);

    ZZ power_of_m(INIT_SIZE, n.size());
    ZZ c_i(INIT_SIZE, n.size());
    ZZ tmp_n = n;

    power(power_of_m, m, d);

    tmp_n = tmp_n - leadingCoef*power_of_m;

    SetCoeff(poly, d, leadingCoef); //The leading coefficient is known...Just find the rest....

    for (long i = d - 1; i >= 0; --i) {
        div(power_of_m, power_of_m, m);
        DivRem(c_i, tmp_n, tmp_n, power_of_m);
        SetCoeff(poly, i, c_i);
    }

    return poly;
}

void PolynomialSelection::tweakRestCoeffsOfKnownNonSkewedNonMonicPolynomial(ZZX &poly, const ZZ& m, bool &changed) {

    RR checkVal = to_RR(m / 2);

    for (long i = 0; i < d - 1; i++) {
        if (poly.rep[i] > to_ZZ(checkVal)) {
            SetCoeff(poly, i, poly.rep[i] - m);
            SetCoeff(poly, i + 1, poly.rep[i + 1] + 1);

            changed = true;
        }
    }
}

bool PolynomialSelection::checkCoeffsSizeToDiscardPolynomial(ZZX &poly, const ZZ& m) {

    RR checkVal = to_RR(m / 2);

    for (long i = 0; i < d - 1; i++) {
        if (abs(poly.rep[i]) > to_ZZ(checkVal)) {
            //cout << "poly.rep[" << i << "] -> " << abs(poly.rep[i]) << " --- m -> " << to_ZZ(checkVal) << endl;
            return true;
        }
    }

    return false;
}

void PolynomialSelection::findBestNonSkewedNonMonicPolynomial(ZZX& f, ZZ& m, const long B, const ZZ& minM, const ZZ& maxM, const double x1, const double x2, const long range, const bool byRating, const long numOfBestRatedPolynomialsToStore, const bool performSieveTest) {
    priorityQueueForRatingPolynomials topRatedPolynomialQueue(ratedPolynomials((bool) byRating));
    //priorityQueueForRatingPolynomials topRatedPolynomialQueue;
    ratedPolynomials polynomial;

    ZZ cofactor;

    vector<long> primesLessThan100;
    getPrimesLessThan100(primesLessThan100);

    //ZZ testM = minM + (maxM - minM) / 2; //We start with an m, which is minM<m<maxM. Let's say in the middle

    double logN = log(n);
    ZZ testM;
    ZZ leadingCoef;
    RR tempPow;

    ZZ aMin;
    ZZ aMax;

    //ZZ aMin = testM / 20000;

    //todo
    /*
        ZZ aMin = testM / 3001; //todo --- Na dw an mporw na to paragw diaforetika (pio dinamika)
        cout << aMin << endl;
        ZZ aMax = testM / 3000;
        aMax = aMin + 10000;
        cout << aMax << endl;
     */


    //todo
    /*
    #define MAX_SHRINKAGE 0.00001
    #define MIN_SHRINKAGE 0.0004
    #define WORST_SHRINKAGE 0.02
     */
    /*
    #define MAX_SHRINKAGE 0.00001
    #define MIN_SHRINKAGE 0.004
    #define WORST_SHRINKAGE 0.02
     */

#ifdef GNFS_DEBUG
    cout << "x1 -> " << x1 << endl;
    cout << "x2 -> " << x2 << endl;
    cout << "range -> " << range << endl;
#endif

    //cout << "exp((d * log(x1) + logN) / (d + 1)) ---> " << exp((d * log(x1) + logN) / (d + 1)) << endl;
    //cout << "exp((d * log(x2) + logN) / (d + 1)) ---> " << exp((d * log(x2) + logN) / (d + 1)) << endl << endl;

    if (range <= 0) {
        aMin = to_ZZ(exp((d * log(x1) + logN) / (d + 1)));
        aMax = to_ZZ(exp((d * log(x2) + logN) / (d + 1)));
    }
    else {
        aMin = to_ZZ(exp((d * log(x1) + logN) / (d + 1)));
        aMax = aMin + range;
    }

#ifdef GNFS_DEBUG
    cout << "aMin ---> " << aMin << endl;
    cout << "aMax ---> " << aMax << endl << endl;
#endif


    //todo 
    //by Jason Papadopoulos...
    //
    // Edw na ftiaksw to ginomeno twn primes gia katallilo kathe fora cofactor.
    // Autos o kwdikas einai apo ton Jason Papadopoulo (msieve)
    if (aMin < to_ZZ(1e8))
        cofactor = 2 * 2 * 3;
    else if (aMin < to_ZZ(1e13))
        cofactor = 2 * 2 * 3 * 3 * 5;
    else if (aMin < to_ZZ(1e15))
        cofactor = 2 * 2 * 2 * 3 * 3 * 3 * 5;
    else if (aMin < to_ZZ(5e17))
        cofactor = 2 * 2 * 2 * 3 * 3 * 3 * 5 * 5 * 7;
    else if (aMin < to_ZZ(1e20))
        cofactor = 2 * 2 * 2 * 3 * 3 * 3 * 5 * 5 * 5 * 7;
    else
        cofactor = 2 * 2 * 2 * 3 * 3 * 3 * 5 * 5 * 5 * 7 * 11;

#ifdef GNFS_DEBUG
    cout << "cofactor ---> " << cofactor << endl;
#endif

    leadingCoef = aMin;

    // The leading coefficient is at least one (1).
    if (leadingCoef == 0)
        leadingCoef++;

    while (leadingCoef <= aMax) {
        /*
         * If the leading coefficient is  not divisible by the coefficient in
         * the disired interval, then update it so that it be divisible.
         */
        while (leadingCoef % cofactor != 0) {
            leadingCoef++; //Now the leading coefficient is divisible by the coefficient
            cout << "Updated leadingCoef -> " << leadingCoef << endl;
        }

        /*
         * Now that cofactor c of a_d divides a_d we can determine m from Eq.(5.3) in Murphy's thesis.
         */

        //todo
        //
        //Na ftiaksw to precision otan kalw RR
        pow(tempPow, to_RR(n) / to_RR(leadingCoef), to_RR(1) / d);
        FloorToZZ(testM, tempPow);


#ifdef GNFS_DEBUG
        /*
         * This cout allows us to see how fast leadingCoef & testM converge.
         * When those two converge, then the iterative process stops.
         */
        cout << "leadingCoef -> " << leadingCoef << "  --- aMax -> " << aMax << "  --- Remaining -> " << aMax - leadingCoef << endl;
#endif

        /*
         * We know what the leading coefficient is, and we want to find the rest of them.
         */
        polynomial.poly = findNonSkewedNonMonicPolynomialWithLeadingCoefKnown(testM, leadingCoef);
        polynomial.m = testM;

        bool changed = false; //Not used in here, but is needed as a parameter in tweakRestCoeffsOfKnownNonSkewedNonMonicPolynomial()
        /*
         * Now we tweak the coefficients. After tweaking all of them are less than m/2.
         */
        tweakRestCoeffsOfKnownNonSkewedNonMonicPolynomial(polynomial.poly, polynomial.m, changed);

        polynomial.alpha = calcAlpha(polynomial.poly, primesLessThan100);

        if (byRating)
            polynomial.rating = ratePolynomial(polynomial.poly, polynomial.alpha, B);
        else
            polynomial.rating = 0; // Now rating=0. Later we will calculate it only for those with the least alpha value.

        // We store only "numOfBestRatedPolynomialsToStore" polynomials in the priority_queue.
        if (topRatedPolynomialQueue.size() < numOfBestRatedPolynomialsToStore) {
            topRatedPolynomialQueue.push(polynomial);
        }
        else {
            /*
             * First we check whether we want to sort by highest rating or least alpha value
             * and then we check whether the new polynomial has higher rating (less alpha) value
             * respectively, than the one on the top of the priority_queue. If it does, then we
             * pop the one on the top and push the new one into the queue. If it doesn't we do not
             * change anything in the queue. 
             */
            if (byRating) {
                if (topRatedPolynomialQueue.top().rating < polynomial.rating) {
                    topRatedPolynomialQueue.pop(); //We pop the polynomial with the worst alpha.
                    topRatedPolynomialQueue.push(polynomial); //and insert the new one with the smaller alpha
                }
            }
            else {
                if (topRatedPolynomialQueue.top().alpha > polynomial.alpha) {
                    topRatedPolynomialQueue.pop(); //We pop the polynomial with the worst alpha.
                    topRatedPolynomialQueue.push(polynomial); //and insert the new one with the smaller alpha
                }
            }
        }
        leadingCoef += cofactor; // Increment the leading coefficient by a step of cofactor
    }

    /*
     * Now from the selected polynomials we have to decide which one we will finally use
     * for the execution of the algorithm. The decision is made either based on the results
     * of some small sieving tests, or on the alpha value of each polynomial.
     */

    /*
     * if performSieveTest == true, then we perform small sieve test on the selected
     * polynomials to find the one that gives the most smooth values.
     *
     * THIS IS NOT IMPLEMENTED YET!!!
     */
    if (performSieveTest) {
        //todo
        //
        // NOT IMPLEMENTED YET !!!
        //
        // Here we will choose the best polynomial from the priority queue
        // through a small sieve test. This is to be implemented in future
        // version.
    }
    else { //else we select the one with the least alpha value, which is the one at the bottom of the priority_queue

        //Anoigw to arxeio gia grapsimo...
        ofstream out(bestPolynomialsFile);

        if (!out.is_open()) {
            cout << "Error opening " << bestPolynomialsFile << ". Exiting..." << endl;
            exit(EXIT_FAILURE);
        }

        ZZX diffF = ZZX::zero();

        while (!topRatedPolynomialQueue.empty()) {


            polynomial = topRatedPolynomialQueue.top();

            diffF = ZZX::zero();
            diffF = diff(polynomial.poly);

            // Here we rate only the polynomials that haven't already been rated.
            if (!byRating && polynomial.rating == 0)
                polynomial.rating = ratePolynomial(polynomial.poly, polynomial.alpha, B);

            //cout << "polynomial -> " << polynomial.poly << " -- m -> " << polynomial.m << " -- alpha -> " << polynomial.alpha << " -- rating -> " << polynomial.rating << endl;

            out << "polynomial -> " << polynomial.poly << endl << "m -> " << polynomial.m << endl << "diffF -> " << diffF << endl << "alpha -> " << polynomial.alpha << endl << "rating -> " << polynomial.rating << endl << endl;



            //We pop the elements of the queue until we reah to the one at the bottom,
            //which is tha one with the least alpha value. This is the polynomials which
            //is finally returned. All the rest are just written to the file.
            if (topRatedPolynomialQueue.size() == 1) {
                //polynomial = topRatedPolynomialQueue.top();

                cout << "---------------------------------------------------------" << endl;
                cout << "---Non-skewed, Non-monic polynomial selection results ---" << endl;
                cout << "---------------------------------------------------------" << endl;
                cout << "-- polynomial -> " << polynomial.poly << endl;
                cout << "-- m -> " << polynomial.m << endl;
                cout << "-- alpha -> " << polynomial.alpha << endl;
                cout << "-- rating -> " << polynomial.rating << endl;
                cout << "---------------------------------------------------------" << endl;

                f = polynomial.poly;
                m = polynomial.m;

                out.close();
            }
            topRatedPolynomialQueue.pop();
        }
    }
}

/*vector<ZZX> PolynomialSelection::getAllPolynomialsZZX(Parameters& params) {

    vector<ZZX> polynomials;
    ZZ minM = params.minM;
    ZZ maxM = params.maxM;
    ZZ m = to_ZZ(0); //Tha krataei to current m gia to opoio tha kalesw tin findCoefficientsForMonicPolynomial()

    for (m = minM; m <= maxM; m++) {
        ZZX coefsForNonTweakedMonic = findMonicPolynomial(m);
        polynomials.push_back(coefsForNonTweakedMonic);
        mValueofEachPolynomial.push_back(m); //Tha einai 1-1 oi times twn polyonimwn kai twn m pou antistoixoun se auta.

        ZZX tweakedPoly = tweakMonicPolynomial(params, coefsForNonTweakedMonic, m);
        polynomials.push_back(tweakedPoly);
        mValueofEachPolynomial.push_back(m); //Tha einai 1-1 oi times twn polyonimwn kai twn m pou antistoixoun se auta.
    }

    return polynomials;
}

ZZX PolynomialSelection::tweakMonicPolynomial(Parameters& params, ZZX& coeffs, ZZ& m) {

    //int n = params.getN();
    int d = params.d;

    //int *tweakedCoeffs = (int *) malloc((int) (d + 1) * sizeof (int));
    ZZX tweakedCoeffs; // = new ZZ[d + 1];
    tweakedCoeffs.rep.SetLength(d + 1);

    //    ZZ *oldCoeffs = coeffs;
    //ZZX oldCoeffs = coeffs; //TO TEST

    //Initialize the tweaked polynomial with the values of the initial polynomial
    for (int k = 0; k < d + 1; k++) {
        tweakedCoeffs.rep[k] = coeffs.rep[k];
    }

    RR checkVal = to_RR(m / 2);
    //cout << "checkVal->" << checkVal << endl;

    for (int i = 0; i < d - 1; i++) {
        if (coeffs.rep[i] > to_ZZ(checkVal)) {
            //cout << "KATI PREPEI NA KANW GIA THN c" << i << endl;
            tweakedCoeffs.rep[i] = tweakedCoeffs.rep[i] - m;
            tweakedCoeffs.rep[i + 1] = tweakedCoeffs.rep[i + 1] + 1;

            //cout << "c" << i << "->" << tweakedCoeffs[i] << endl;
            //cout << "c" << i + 1 << "->" << tweakedCoeffs[i + 1] << endl;
        }
    }

    return tweakedCoeffs;
}

ZZX PolynomialSelection::testAllPolynomialsAndReturnTheBestOne(vector<ZZX>& allPolynomials) {

    //Pros to paron apla tha epistrefei ena polyonymo, estw to prwto.

    int best = 34; //Se poia thesi vrisketai to kalytero polyonimo (best=35 || best=34 => m=31)

    ZZX polynomialToReturn = allPolynomials.at(best);

    //cout << "m=" << mValueofEachPolynomial.at(3) << endl;
    m = mValueofEachPolynomial.at(best);

    return polynomialToReturn;
}
 */

bool PolynomialSelection::findBestPolynomial(ZZX& f, ZZ& m, const ZZ& minM, const ZZ& maxM, const long B, const bool byRating, const long numOfBestRatedPolynomialsToStore, const bool performSieveTest) {
    vector<ratedPolynomials> topRatedPolynomialsVector;
    topRatedPolynomialsVector.reserve(numOfBestRatedPolynomialsToStore);

    //cout << endl << endl << "numOfPolynomials -> " << numOfPolynomials << " && maxQueueSize -> " << maxQueueSize << endl << endl;

    //We now find the best rated polynomials
    findTopRatedPolynomials(topRatedPolynomialsVector, minM, maxM, B, byRating, numOfBestRatedPolynomialsToStore);

    if (topRatedPolynomialsVector.empty()) {
        return false;
    }

    ratedPolynomials bestRatedPolynomial;

    /*
     * Now from the selected polynomials we have to decide which one we will finally use
     * for the execution of the algorithm. The decision is made either based on the results
     * of some small sieving tests, or on the alpha value of each polynomial.
     */

    /*
     * if performSieveTest == true, then we perform small sieve test on the selected
     * polynomials to find the one that gives the most smooth values.
     *
     * THIS IS NOT IMPLEMENTED YET!!!
     */
    if (performSieveTest) {
        //todo
        //
        // NOT IMPLEMENTED YET !!!
        //
        // Here we will choose the best polynomial from the vector
        // through a small sieve test. This is to be implemented in future
        // version.
    }
    else { //else we select the one with the highest rating value, which is the one at the begining of the vector

        /*
         * Now we store in a file all the best rated polynomials found by the process.
         * Finally the one with the highest rating is selected as the one that will
         * be used for the algorithm. The rest in the file are stored for manual use,
         * in case the best rated one fails to give us the number of relations needed.
         */
        ofstream out(bestPolynomialsFile);

        if (!out.is_open()) {
            cout << "Error opening " << bestPolynomialsFile << ". Exiting..." << endl;
            exit(EXIT_FAILURE);
        }

        ZZX diffF = ZZX::zero();

        vector<ratedPolynomials>::iterator it;
        for (it = topRatedPolynomialsVector.begin(); it != topRatedPolynomialsVector.end(); ++it) {

            bestRatedPolynomial = *it;

            diffF = ZZX::zero();
            diffF = diff(bestRatedPolynomial.poly);


            //todo
            //DE XREIAZETAI NA KANW PALI RATING GIATI HDH EXW KANEI GIA OLA TA POLYWNUMA POU THA ERTHOUN EDW.
            //bestRatedPolynomial.rating = ratePolynomial(bestRatedPolynomial.poly, bestRatedPolynomial.alpha, B);

            out << "polynomial -> " << bestRatedPolynomial.poly << endl << "m -> " << bestRatedPolynomial.m << endl << "diffF -> " << diffF << endl << "alpha -> " << bestRatedPolynomial.alpha << endl << "rating -> " << bestRatedPolynomial.rating << endl << endl;

        }

        out.close(); //We now close the file



        //The best rated polynomial is the one in the begining of the vector.
        bestRatedPolynomial = topRatedPolynomialsVector[0];

        cout << "polynomial -> " << bestRatedPolynomial.poly << " -- m -> " << bestRatedPolynomial.m << " -- alpha -> " << bestRatedPolynomial.alpha << " -- rating -> " << bestRatedPolynomial.rating << endl;

        f = bestRatedPolynomial.poly;
        m = bestRatedPolynomial.m;
    }

    return true;
}

void PolynomialSelection::findTopRatedPolynomials(vector<ratedPolynomials>& topRatedPolynomialsVector, const ZZ& minM, const ZZ& maxM, const long B, const bool byRating, long numOfBestRatedPolynomialsToStore) {
    ZZ m; // This var stores the temp value of m for which we will call findMonicPolynomial()
    //priorityQueueForRatingPolynomials topRatedPolynomialQueue;
    priorityQueueForRatingPolynomials topRatedPolynomialQueue(ratedPolynomials((bool) byRating));

    ratedPolynomials polynomial;

    long maxNumOfPolynomials = 100; //The number of polynomials that will be tested

    //We store the vector with the primes p, with p<100, to use it in alpha's calculation.
    vector<long> primesLessThan100;
    getPrimesLessThan100(primesLessThan100);

    /*
     * If the diferrence of "maxM - minM" is less than the number of the polynomials
     * we want to store, then we (obviously) want less!!! ;)
     */
    if ((maxM - minM) < maxNumOfPolynomials)
        maxNumOfPolynomials = to_long(maxM - minM);
    /*
         if ((maxM - minM) < numOfBestRatedPolynomialsToStore)
         numOfBestRatedPolynomialsToStore = to_long(maxM - minM);
     */



    /*
     * We start the search process for m = maxM, which is the m-value returned
     * from the base-m expansion. M is reduced by one each time, until we have
     * found "numOfBestRatedPolynomialsToStore". The main idea is that the best
     * polynomial will be near and below that initial m.
     */
    m = maxM;

    bool changed = false;

    //for (long i = 0; i < numOfBestRatedPolynomialsToStore; i++) {
    for (long i = 0; i < maxNumOfPolynomials; i++) {
        polynomial.poly = findMonicPolynomial(m);
        polynomial.m = m;
        polynomial.alpha = calcAlpha(polynomial.poly, primesLessThan100);

        for (long tweakIDX = 0; tweakIDX < 2; tweakIDX++) {

            changed = false;

            /*
             * If we are in the second iteration, then we tweak the last selected polynomial.
             */
            if (tweakIDX == 1) {
                tweakRestCoeffsOfKnownNonSkewedNonMonicPolynomial(polynomial.poly, polynomial.m, changed);

                //cout << "changed -> " << changed << endl;

                /*
                 * If the tweaking process didn't change any coefficient then do
                 * not put it again into the priority_queu, just break.
                 */
                if (!changed)
                    break;
            }



            //if true then the polynomial is not good (a_i > m/2), so do not rate it.
            if (checkCoeffsSizeToDiscardPolynomial(polynomial.poly, polynomial.m)) {
                //cout << "NOT GOOD POLYNOMIAL -> " << polynomial.poly << "  with m: " << polynomial.m << endl;
                continue;
            }
            else {
                if (byRating)
                    polynomial.rating = ratePolynomial(polynomial.poly, polynomial.alpha, B);
                else
                    polynomial.rating = 0; // Now rating=0. Later we will calculate it only for those with the least alpha value.

                /*
                 * Now we check the rating or the alpha value of each polynomial and we store the
                 * "numOfBestRatedPolynomialsToStore" best (lower alpha value or higher rating value) in a
                 * priority_queue.
                 */
                if (topRatedPolynomialQueue.size() < numOfBestRatedPolynomialsToStore) {
                    topRatedPolynomialQueue.push(polynomial);
                }
                else {
                    /*
                     * First we check whether we want to sort by highest rating or least alpha value
                     * and then we check whether the new polynomial has higher rating (less alpha) value
                     * respectively, than the one on the top of the priority_queue. If it does, then we
                     * pop the one on the top and push the new one into the queue. If it doesn't we do not
                     * change anything in the queue.
                     */
                    if (byRating) {
                        if (topRatedPolynomialQueue.top().rating < polynomial.rating) {
                            topRatedPolynomialQueue.pop(); //We pop the polynomial with the worst alpha.
                            topRatedPolynomialQueue.push(polynomial); //and insert the new one with the smaller alpha
                        }
                    }
                    else {
                        if (topRatedPolynomialQueue.top().alpha > polynomial.alpha) {
                            topRatedPolynomialQueue.pop(); //We pop the polynomial with the worst alpha.
                            topRatedPolynomialQueue.push(polynomial); //and insert the new one with the smaller alpha
                        }
                    }
                }
            }

        }
        m--;
    }




    /*
     * We create the polynomials for each m with minM < m < maxM and we store in a
     * priority_queue those with the least alpha value. The number of the stored
     * polynomials is equal to maxQueueSize.
     */
    /*
         for (m = minM; m <= maxM; m++) {
         polynomial.poly = findMonicPolynomial(m);
         polynomial.m = m;
         polynomial.alpha = calcAlpha(polynomial.poly, primesLessThan100);
         polynomial.rating = 0; // Now rating=0. Later we will calculate it only for those with the least alpha value.
	 
         if (topRatedPolynomialQueue.size() < numOfBestRatedPolynomialsToStore) {
         topRatedPolynomialQueue.push(polynomial);
         }
         else {
         // We check whether the new polynomial has smaller alpha value than the
         // one on the top of the priority_queue. If it does, then we pop the one
         //  on the top and push the new one into the queue. If it doesn't we do not
         // change anything in the queue.
         if (topRatedPolynomialQueue.top().alpha > polynomial.alpha) {
         topRatedPolynomialQueue.pop(); //We pop the polynomial with the worst alpha.
         topRatedPolynomialQueue.push(polynomial); //and insert the new one with the smaller alpha
         }
         }
         }
     */

    /*
     * Now we rate each polynomial and store it in the vector which will
     * finally be sorted by each polynomial's rating.
     */
    while (!topRatedPolynomialQueue.empty()) {
        //cout << topRatedPolynomialQueue.top().poly <<" - " << topRatedPolynomialQueue.top().m<< " - " << topRatedPolynomialQueue.top().alpha<<endl;

        polynomial = topRatedPolynomialQueue.top();

        //Now we can rate each polynomial (if we haven't already rated it).
        if (!byRating)
            polynomial.rating = ratePolynomial(polynomial.poly, polynomial.alpha, B);

        //and store it in the vector
        topRatedPolynomialsVector.push_back(polynomial);

        topRatedPolynomialQueue.pop();
    }

    //TO DO
    /*
         cout << endl << endl;
         vector<ratedPolynomials>::iterator it;
         cout << "myvector contains before sorting:" << endl;
         for (it = topRatedPolynomialsVector.begin(); it != topRatedPolynomialsVector.end(); ++it) {
         polynomial = *it;
         cout << "polynomial.poly -> " << polynomial.poly << " ----- alpha -> " << polynomial.alpha << " ----- rating -> " << polynomial.rating << endl;
         }
     */


    /*
     * We now sort the polynomials into the vector based on their rating. The one with the highest rating
     * is stored at the begining of the vector.
     */
    sort(topRatedPolynomialsVector.begin(), topRatedPolynomialsVector.end(), checkRating);


    //TO DO
    /*
         cout << endl << endl;
         //vector<ratedPolynomials>::iterator it;
         cout << "myvector contains after sorting:" << endl;
         for (it = topRatedPolynomialsVector.begin(); it != topRatedPolynomialsVector.end(); ++it) {
         polynomial = *it;
         //cout << " " << *it;
         cout << "polynomial.poly -> " << polynomial.poly << " ----- alpha -> " << polynomial.alpha << " ----- rating -> " << polynomial.rating << endl;
         }
         cout << endl << endl;
     */
}

bool PolynomialSelection::findBestPolynomial_RR(ZZX& f, ZZ& m, const ZZ& minM, const ZZ& maxM, const long B, const long numOfBestRatedPolynomialsToStore, const bool performSieveTest) {
    vector<ratedPolynomials_RR> topRatedPolynomialsVector;
    topRatedPolynomialsVector.reserve(numOfBestRatedPolynomialsToStore);

    //cout << endl << endl << "numOfPolynomials -> " << numOfPolynomials << " && maxQueueSize -> " << maxQueueSize << endl << endl;

    //We now find the best rated polynomials
    findTopRatedPolynomials_RR(topRatedPolynomialsVector, minM, maxM, B, numOfBestRatedPolynomialsToStore);

    if (topRatedPolynomialsVector.empty()) {
        return false;
    }

    ratedPolynomials_RR bestRatedPolynomial;

    /*
     * Now from the selected polynomials we have to decide which one we will finally use
     * for the execution of the algorithm. The decision is made either based on the results
     * of some small sieving tests, or on the alpha value of each polynomial.
     */

    /*
     * if performSieveTest == true, then we perform small sieve test on the selected
     * polynomials to find the one that gives the most smooth values.
     *
     * THIS IS NOT IMPLEMENTED YET!!!
     */
    if (performSieveTest) {
        //todo
        //
        // NOT IMPLEMENTED YET !!!
        //
        // Here we will choose the best polynomial from the vector
        // through a small sieve test. This is to be implemented in future
        // version.
    }
    else { //else we select the one with the highest rating value, which is the one at the begining of the vector
        //The best rated polynomial is the one in the begining of the vector.
        bestRatedPolynomial = topRatedPolynomialsVector[0];

        cout << "polynomial -> " << bestRatedPolynomial.poly << " -- m -> " << bestRatedPolynomial.m << " -- alpha -> " << bestRatedPolynomial.alpha << " -- rating -> " << bestRatedPolynomial.rating << endl;

        f = bestRatedPolynomial.poly;
        m = bestRatedPolynomial.m;
    }

    return true;
}

void PolynomialSelection::findTopRatedPolynomials_RR(vector<ratedPolynomials_RR>& topRatedPolynomialsVector, const ZZ& minM, const ZZ& maxM, const long B, const long numOfBestRatedPolynomialsToStore) {
    ZZ m; // This var stores the temp value of m for which we will call findMonicPolynomial()
    priorityQueueForRatingPolynomials_RR topRatedPolynomialQueue;

    ratedPolynomials_RR polynomial;

    //We store the vector with the primes p, with p<100, to use it in alpha's calculation.
    vector<long> primesLessThan100;
    getPrimesLessThan100(primesLessThan100);

    /*
     * We create the polynomials for each m with minM < m < maxM and we store in a
     * priority_queue those with the least alpha value. The number of the stored
     * polynomials is equal to maxQueueSize.
     */
    for (m = minM; m <= maxM; m++) {
        polynomial.poly = findMonicPolynomial(m);
        polynomial.m = m;
        polynomial.alpha = calcAlpha_RR(polynomial.poly, primesLessThan100);
        polynomial.rating = to_RR(0); // Now rating=0. Later we will calculate it only for those with the least alpha value.

        if (topRatedPolynomialQueue.size() < numOfBestRatedPolynomialsToStore) {
            topRatedPolynomialQueue.push(polynomial);
        }
        else {
            /*
             * We check whether the new polynomial has smaller alpha value than the
             * one on the top of the priority_queue. If it does, then we pop the one
             * on the top and push the new one into the queue. If it doesn't we do not
             * change anything in the queue.
             */
            if (topRatedPolynomialQueue.top().alpha > polynomial.alpha) {
                topRatedPolynomialQueue.pop(); //We pop the polynomial with the worst alpha.
                topRatedPolynomialQueue.push(polynomial); //and insert the new one with the smaller alpha
            }
        }
    }

    /*
     * Now we rate each polynomial and store it in the vector which will
     * finally be sorted by each polynomial's rating.
     */
    while (!topRatedPolynomialQueue.empty()) {
        //cout << topRatedPolynomialQueue.top().poly <<" - " << topRatedPolynomialQueue.top().m<< " - " << topRatedPolynomialQueue.top().alpha<<endl;

        polynomial = topRatedPolynomialQueue.top();

        //Now we can rate each polynomial
        polynomial.rating = ratePolynomial_RR(polynomial.poly, polynomial.alpha, B);

        //and store it in the vector
        topRatedPolynomialsVector.push_back(polynomial);

        topRatedPolynomialQueue.pop();
    }

    //TO DO
    /*
            cout << endl << endl;
            vector<ratedPolynomials>::iterator it;
            cout << "myvector contains before sorting:" << endl;
            for (it = topRatedPolynomialsVector.begin(); it != topRatedPolynomialsVector.end(); ++it) {
                polynomial = *it;
                cout << "polynomial.poly -> " << polynomial.poly << " ----- alpha -> " << polynomial.alpha << " ----- rating -> " << polynomial.rating << endl;
            }
     */


    /*
     * We now sort the polynomials into the vector based on their rating. The one with the highest rating
     * is stored at the begining of the vector.
     */
    sort(topRatedPolynomialsVector.begin(), topRatedPolynomialsVector.end(), checkRating_RR);


    //TO DO
    /*
            cout << endl << endl;
            //vector<ratedPolynomials>::iterator it;
            cout << "myvector contains after sorting:" << endl;
            for (it = topRatedPolynomialsVector.begin(); it != topRatedPolynomialsVector.end(); ++it) {
                polynomial = *it;
                //cout << " " << *it;
                cout << "polynomial.poly -> " << polynomial.poly << " ----- alpha -> " << polynomial.alpha << " ----- rating -> " << polynomial.rating << endl;
            }
            cout << endl << endl;
     */
}

/*
void PolynomialSelection::storeAllPolynomialsInFile(vector<ZZX>& allPolynomials, vector<ZZ>& mValueofEachPolynomial) {

    ofstream output;
    output.open("polynomials.txt", ios::out); // open the streams

    output << "#Polynomial" << " \t" << "m" << "\t Tweaked?" << endl << endl;

    int cnt = 0;
    vector<ZZ>::iterator it2 = mValueofEachPolynomial.begin();

    for (vector<ZZX>::iterator it = allPolynomials.begin(); it != allPolynomials.end(); ++it, it2++) {

        string tweaked = "no";
        if (cnt % 2 == 0)
            tweaked = "no";
        else
            tweaked = "yes";


        output << *it << " \t" << *it2 << "\t" << tweaked << endl;
        cnt++;
    }

    output << endl << "Number of Polynomials->" << cnt << endl;

    output.close();
}
 */

void PolynomialSelection::getPrimesLessThan100(vector<long>& primesLessThan100) {
    PrimeSeq s; // a class for quickly generating primes in sequence
    long p;

    p = s.next(); //get the first prime

    while (p < 100) {
        primesLessThan100.push_back(p);
        p = s.next();
    }
}

/*
 * Version with more precision. Calculates a.
 */
RR PolynomialSelection::calcAlpha_RR(ZZX &f, vector<long> &primes) {
    ZZ disc = discriminant(f);
    RR cont_p = to_RR(0);
    RR alpha = to_RR(0);

    //This is the sample bound for the poorly-behaved primes (p|Discriminant)
    //and is used to evaluate Eq. (3.8) in Murphy's thesis
    int bound = 10000;
    RR p;

    //We compute the contribution of each prime p
    for (vector<long>::iterator iter = primes.begin(); iter != primes.end(); ++iter) {
        p = to_RR(*iter);

        cont_p = computeContP_RR(f, p, disc, bound);

        //Primes p at which f has no roots are omitted [Murphy $3.2.2]
        //if (cont_p > 0)
        alpha += (inv(p - 1) - cont_p) * log(p); //Quantify root properties (a) [Murphy $3.2.3]
    }
    return alpha;
}

double PolynomialSelection::calcAlpha(ZZX &f, vector<long> &primes) {
    ZZ disc = discriminant(f);
    double cont_p = 0;
    double alpha = 0;

    //This is the sample bound for the poorly-behaved primes (p|Discriminant)
    //and is used to evaluate Eq. (3.8) in Murphy's thesis
    int bound = 10000;
    long p;

    //We compute the contribution of each prime p
    for (vector<long>::iterator iter = primes.begin(); iter != primes.end(); ++iter) {
        p = *iter;

        cont_p = computeContP(f, p, disc, bound);

        //Primes p at which f has no roots are omitted [Murphy $3.2.2]
        //if (cont_p > 0)
        alpha += (1.0 / (p - 1) - cont_p) * log(p); //Quantify root properties (a) [Murphy $3.2.3]
        //alpha += (inv(p - 1) - cont_p) * log(p); //Quantify root properties (a) [Murphy $3.2.3]
    }

    return alpha;
}

/*
 * Version with better precision. Computes the contribution of each prime.
 */
RR PolynomialSelection::computeContP_RR(ZZX& f, RR& prime, ZZ& disc, int bound) {
    ZZ p = to_ZZ(prime);
    ZZ a, b, value;
    ZZ a_d = LeadCoeff(f); //The leading coefficient, which is used to find the projective roots
    ZZ upperBound = to_ZZ(pow(to_RR(bound), to_RR(2)));

    long q_p = 0;
    value = to_ZZ(0);
    RR cont_p = to_RR(0);

    //If p is a poorly-behaved prime (p|Discriminant), then we estimate p's
    //contribution using Eq (3.8) in Murphy's thesis $3.2.1.
    if (disc % p == 0) {
        for (long i = 1; i < bound; i++) {
            do {
                a = RandomBnd(upperBound);
                b = RandomBnd(upperBound);
            } while (GCD(a, b) != 1 || a == 0 || b == 0);


            evalHomogeneousF(value, f, a, b);

            //evalHomogeneousF(value, f, a, b);
            while (value % p == 0) {
                value = value / p;
                cont_p++;
            }
        }
        cont_p = cont_p / bound;
    }
    else { //else if p is a well-behaved prime we compute its full contribution
        zz_p::init(to_long(p));
        zz_pX fx = to_zz_pX(f); //fmod(p)

        // If the polynomial is not monic, then convert it into its monic associate
        if (a_d > 1) {
            /*
             * If the polynomial is non-monic, then convert it into monic and then
             * find its roots modulo p.
             *
             * LeadCoeff(fx) must be invertible modulo p, which is, here, always the case.
             * This happens because the polynomial fx is reduced modulo p, which means
             * that all the coefficients are less than p. In our case, p is prime,
             * and this means that p and LeadCoeff(fx) are relatively prime
             * [ GCD( p, LeadCoeff(fx) ) = 1 ], which is the case for LeadCoeff(fx)
             * to be invertible.
             */
            MakeMonic(fx);
        }

#ifdef GNFS_FACTOR_BASES_USE_FINDROOTS

        zz_pX X;
        SetX(X);
        zz_pX g = PowerXMod(p, fx) - X;

        vec_zz_p roots;
        FindRoots(roots, GCD(fx, g));

        long length = roots.length();
        for (long i = 0; i < length; ++i) {
            q_p++;
        }
#else
        vec_pair_zz_pX_long factors;
        CanZass(factors, fx);

        //We count the real roots of fmod(p)
        for (int i = 0; i < factors.length(); i++) {
            if (deg(factors[i].a) == 1) {
                q_p++;
            }
        }
#endif

        //Now we check for contribution from Projective Roots
        if (a_d % p == 0)
            q_p++;

        //Primes p at which f has no roots are omitted [Murphy $3.2.2]
        //if (q_p > 0)
        //cont_p = to_RR(q_p) / (to_RR(p) - to_RR(1)); // Murphy $3.2.2, Eq.(3.10)
        cont_p = (to_RR(q_p) * to_RR(p)) / (pow(to_RR(p), to_RR(2)) - to_RR(1)); // Murphy $3.2.2, Eq.(3.11)
    }
    return cont_p;
}

/*
 * Version with less precision (double). Computes the contribution of each prime.
 */
double PolynomialSelection::computeContP(ZZX& f, long p, ZZ& disc, int bound) {
    ZZ a, b, value;
    ZZ a_d = LeadCoeff(f); //The leading coefficient, which is used to find the projective roots
    ZZ upperBound = to_ZZ(pow(to_RR(bound), to_RR(2)));

    long q_p = 0;
    value = to_ZZ(0);
    double cont_p = 0;


    //If p is a poorly-behaved prime (p|Discriminant), then we estimate p's
    //contribution using Eq (3.8) in Murphy's thesis $3.2.1.
    if (disc % p == 0) {
        for (long i = 1; i < bound; i++) {
            do {
                a = RandomBnd(upperBound);
                b = RandomBnd(upperBound);
            } while (GCD(a, b) != 1 || a == 0 || b == 0);

            evalHomogeneousF(value, f, a, b);
            //evalHomogeneousF(value, f, a, b);
            while (value % p == 0) {
                value = value / p;
                cont_p++;
            }
        }
        cont_p = cont_p / bound;

    }
    else { //else if p is a well-behaved prime we compute its full contribution
        zz_p::init(p);
        zz_pX fx = to_zz_pX(f); //fmod(p)

        //vec_pair_zz_pX_long factors;
        //CanZass(factors, fx);

        // If the polynomial is not monic, then convert it into its monic associate
        if (a_d > 1) {
            /*
             * If the polynomial is non-monic, then convert it into monic and then 
             * find its roots modulo p.
             *
             * LeadCoeff(fx) must be invertible modulo p, which is, here, always the case.
             * This happens because the polynomial fx is reduced modulo p, which means
             * that all the coefficients are less than p. In our case, p is prime,
             * and this means that p and LeadCoeff(fx) are relatively prime
             * [ GCD( p, LeadCoeff(fx) ) = 1 ], which is the case for LeadCoeff(fx)
             * to be invertible.
             */
            MakeMonic(fx);
        }

#ifdef GNFS_FACTOR_BASES_USE_FINDROOTS

        zz_pX X;
        SetX(X);
        zz_pX g = PowerXMod(p, fx) - X;

        vec_zz_p roots;
        FindRoots(roots, GCD(fx, g));

        long length = roots.length();
        for (long i = 0; i < length; ++i) {
            q_p++;
        }
#else
        vec_pair_zz_pX_long factors;
        CanZass(factors, fx);

        //We count the real roots of fmod(p)
        for (int i = 0; i < factors.length(); i++) {
            if (deg(factors[i].a) == 1) {
                q_p++;
            }
        }
#endif

        //Now we check for contribution from Projective Roots
        if (a_d % p == 0)
            q_p++;

        //Primes p at which f has no roots are omitted [Murphy $3.2.2]
        //if (q_p > 0)
        //cont_p = q_p / (p - 1); // Murphy $3.2.2, Eq.(3.10)
        cont_p = ((q_p * p) / (pow(p, 2.0) - 1.0)); // Murphy $3.2.2, Eq.(3.11)

        //cont_p = to_double((to_RR(q_p) * to_RR(p)) / (pow(to_RR(p), to_RR(2)) - to_RR(1))); // Murphy $3.2.2, Eq.(3.11)
        //cont_p = to_RR(q_p) / (to_RR(p) - to_RR(1)); // Murphy $3.2.2, Eq.(3.10)
        //cont_p = q_p / (p - 1); // Murphy $3.2.2, Eq.(3.10)
    }
    return cont_p;
}

/*
 * Compute the rate of a polynomial, based on Murphy's thesis (eq.5.6 p.86)
 */
double PolynomialSelection::ratePolynomial(ZZX& f, double alpha, long B) {
    double rating = 0;
    long K = 1000;

    double u_f = 0;
    double theta = 0;
    double pi = M_PI;
    double logB = log(B);
    long i = 0;
    double COS = 0;
    double SIN = 0;
    double res = 0;

    for (i = 0; i < K; i++) {
        theta = (pi / K) * (i - 0.5);

        COS = cos(theta);
        SIN = sin(theta);

        res = evalHomogeneousF(f, COS, SIN); //todo
        //res = evalHomogeneousF(f, COS, SIN); //todo

        u_f = (log(fabs(res)) + alpha) / (logB + d);

        rating += dickmanRho(u_f);
    }

    return rating / K;
}

/*
 * Compute the rate of a polynomial, based on Murphy's thesis (eq.5.6 p.86)
 */
RR PolynomialSelection::ratePolynomial_RR(ZZX &f, RR &alpha, long B) {
    RR rating = to_RR(0);

    RR u_f = to_RR(0);
    RR theta = to_RR(0);
    RR pi = ComputePi_RR();
    RR logB = log(to_RR(B));
    long K = 1000;

    for (long i = 0; i < K; i++) {
        theta = (pi / K) * (i - 0.5);

        evalHomogeneousF(u_f, f, cos(theta), sin(theta));
        //evalHomogeneousF(u_f, f, cos(theta), sin(theta));//todo
        u_f = (log(fabs(u_f)) + alpha) / (logB + d);

        rating += dickmanRho_RR(u_f);
    }


    return rating / K;
}

double PolynomialSelection::dickmanRho(double& u) {
    /*
     * rho = 1, when u <= 1, as in Back & Peralta, Eq. (2.1).
     */
    if (u < 1.0)
        return 1.0;

    if (u > 50)
        return 0.0;

    int i = 0;
    int j = 0;
    int tmp_k = 0;
    int numOfCoefs = 55;
    int k = ceil(u);
    //double zeta = tmp_k - u;
    double zeta = k - u;
    double *c = new double[numOfCoefs];
    double *tmp_c = new double[numOfCoefs];
    double rho = 0;
    double tmpZeta = 1;

    /*
     * Initialization of the coefficients' matrix
     */
    if (k < 2) {//If k<2, then we have c[0]=1, c[i]=0, for i>=1 [Back & Peralta, Eq. (4.1)]
        c[0] = 1;

        for (i = 1; i < numOfCoefs; i++) {
            c[i] = 0;
        }
        tmp_k = 1;
    }
    else { //else if k>=2 c[0]=1-log2, c[i]=1/(i*2^i), for i>=1 [Back & Peralta, Eq. (4.1)]
        c[0] = 1 - log(2);

        for (i = 1; i < numOfCoefs; i++) {
            c[i] = 1 / (i * pow(2.0, i));
        }
        tmp_k = 2;
    }

    while (tmp_k < k) {
        for (i = 1; i < numOfCoefs; i++) {
            tmp_c[i] = 0; //Always initialize the tmp array

            for (j = 0; j < i; j++) {
                tmp_c[i] += c[j] / (i * pow(tmp_k + 1.0, i - j)); //[Back & Peralta, Eq. (4.2)]
            }
        }

        for (j = 1; j < numOfCoefs; j++)
            tmp_c[0] += tmp_c[j] / (j + 1);

        tmp_c[0] /= tmp_k; //[Back & Peralta, Eq. (4.3)]

        /*
         * Now update the array with the final results
         */
        for (i = 0; i < numOfCoefs; i++)
            c[i] = tmp_c[i];

        tmp_k++;
    }


    /*
     * Now we can compute the value of rho(u) using the (possibly) updated array c
     */
    for (i = 0; i < numOfCoefs; i++) {
        //rho += c[i]*(power(zeta, i));
        rho += c[i] * tmpZeta;
        tmpZeta *= zeta;
    }

    return rho;
}

RR PolynomialSelection::dickmanRho_RR(RR &u) {
    /*
     * rho = 1, when u <= 1, as in Back & Peralta, Eq. (2.1).
     */
    if (u < to_RR(1)) {
        return to_RR(1.0);
    }

    if (u > to_RR(50))
        return to_RR(0.0);

    long i = 0;
    long j = 0;
    long tmp_k = 0;
    long numOfCoefs = 55;
    RR k = ceil(u);
    //RR zeta = tmp_k - u;
    RR zeta = k - u;
    RR *c = new RR[numOfCoefs];
    RR *tmp_c = new RR[numOfCoefs];
    RR rho = to_RR(0);
    RR tmpZeta = to_RR(1);

    RR tmpPow;


    /*
     * Initialization of the coefficients' matrix
     */
    if (k < 2) {//If k<2, then we have c[0]=1, c[i]=0, for i>=1 [Back & Peralta, Eq. (4.1)]
        c[0] = to_RR(1);

        for (i = 1; i < numOfCoefs; i++) {
            c[i] = to_RR(0);
        }
        tmp_k = 1;
    }
    else { //else if k>=2 c[0]=1-log2, c[i]=1/(i*2^i), for i>=1 [Back & Peralta, Eq. (4.1)]
        c[0] = to_RR(1) - log(to_RR(2));

        for (i = 1; i < numOfCoefs; i++) {
            c[i] = inv(to_RR(i) * power(to_RR(2), i));
        }
        tmp_k = 2;
    }

    while (tmp_k < k) {
        for (i = 1; i < numOfCoefs; i++) {

            tmp_c[i] = to_RR(0); //Always initialize the tmp array

            tmpPow = power(to_RR(tmp_k + 1), i);

            for (j = 0; j < i; j++) {
                tmp_c[i] += c[j] / (i * tmpPow); //[Back & Peralta, Eq. (4.2)]
                tmpPow /= to_RR(tmp_k + 1);
            }
        }

        for (j = 1; j < numOfCoefs; j++)
            tmp_c[0] += tmp_c[j] / to_RR(j + 1);

        tmp_c[0] /= tmp_k; //[Back & Peralta, Eq. (4.3)]

        /*
         * Now update the array with the final results
         */
        for (i = 0; i < numOfCoefs; i++)
            c[i] = tmp_c[i];

        tmp_k++;
    }

    /*
     * Now we can compute the value of rho(u) using the (possibly) updated array c
     */
    for (i = 0; i < numOfCoefs; i++) {
        //rho += c[i]*(power(zeta, i));
        rho += c[i] * tmpZeta;
        tmpZeta *= zeta;
    }

    return rho;
}
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


#ifndef TIMER_H_INCLUDED
#define TIMER_H_INCLUDED
#include <exception>

class timer {
public:

    class timer_exception : public std::exception {
        virtual const char *what() const throw ();
    };
    timer();
    ~timer();
    double time() const;
    void reset();
    bool suspend();
    bool resume();
    double elapsed();
private:
    bool suspended_;
    double last_tick_;
    double suspend_time_;
    double base_time_;
    double freq_;
};

#endif


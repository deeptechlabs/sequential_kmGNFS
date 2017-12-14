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


#ifndef _GNFS_RELATION_H
#define	_GNFS_RELATION_H

#include <NTL/vec_long.h>
#include <NTL/vector.h>

#include <GNFS/gnfs_global.h>
#include <pair_intType_intType.h>
#include <triple_intType.h>

NTL_CLIENT

class Relation {
public:

    pair_intType_intType relPair; //(a,b)
    bool sign; //true = negative, false = positive

    vec_pair_intType_intType rationalFactors; //each pair consists of a prime factor and its exponent that occurs in the factorization of F1(a,b)
    vec_triple_intType algebraicFactors; //each triple consists of a prime factor, the corresponding r (p,r) and its exponent that occurs in the factorization of F2(a,b)

    vec_long rationalFactorsIndex;
    vec_long algebraicFactorsIndex;

public:

    void softReset() {

        rationalFactors.SetLength(0L);
        algebraicFactors.SetLength(0L);
        rationalFactorsIndex.SetLength(0L);
        algebraicFactorsIndex.SetLength(0L);
    }

    void appendRationalFactor(const pair_intType_intType& pair, const long factorBaseIndex) {

        append(rationalFactors, pair);
        append(rationalFactorsIndex, factorBaseIndex);
    }

    void appendAlgebraicFactor(const triple_intType& triple, const long factorBaseIndex) {

        append(algebraicFactors, triple);
        append(algebraicFactorsIndex, factorBaseIndex);
    }


private:

};

inline long operator==(const Relation& a, const Relation& b) {
    return a.relPair == b.relPair; //we only check if pair (a,b) is the same and we ignore other info such as sign, factors etc
}

inline long operator!=(const Relation& a, const Relation& b) {
    return !(a == b);
}

std::istream & operator>>(std::istream&, Relation&);
std::ostream & operator<<(std::ostream&, const Relation&);


//vec_Relation type declaration
NTL_vector_decl(Relation, vec_Relation)

NTL_io_vector_decl(Relation, vec_Relation)

NTL_eq_vector_decl(Relation, vec_Relation)

//vec_vec_Relation type declaration
NTL_vector_decl(vec_Relation, vec_vec_Relation)

NTL_eq_vector_decl(vec_Relation, vec_vec_Relation)

NTL_io_vector_decl(vec_Relation, vec_vec_Relation)

#endif	/* _GNFS_RELATION_H */


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


#include <GNFS/Relation.h>

NTL_CLIENT

std::istream & operator>>(std::istream& in, Relation& rel) {

    long c;

    if (!in) NTL_NNS Error("bad Relation object input");

    c = in.peek();
    while (NTL_NNS IsWhiteSpace(c)) {
        in.get();
        c = in.peek();
    }

    if (c != '[')
        NTL_NNS Error("bad Relation object input");

    in.get();

    if (!(in >> rel.sign))
        NTL_NNS Error("bad Relation object input");
    if (!(in >> rel.relPair))
        NTL_NNS Error("bad Relation object input");
    if (!(in >> rel.rationalFactors))
        NTL_NNS Error("bad Relation object input");
    if (!(in >> rel.algebraicFactors))
        NTL_NNS Error("bad Relation object input");
    if (!(in >> rel.rationalFactorsIndex))
        NTL_NNS Error("bad Relation object input");
    if (!(in >> rel.algebraicFactorsIndex))
        NTL_NNS Error("bad Relation object input");

    c = in.peek();
    while (NTL_NNS IsWhiteSpace(c)) {
        in.get();
        c = in.peek();
    }

    if (c != ']')
        NTL_NNS Error("bad Relation object input");

    in.get();

    return in;
}

std::ostream & operator<<(std::ostream& out, const Relation& rel) {

    return out << '[' << rel.sign << rel.relPair << ' ' << rel.rationalFactors << ' ' << rel.algebraicFactors << ' ' << rel.rationalFactorsIndex << ' ' << rel.algebraicFactorsIndex << ']';
}


//vec_Relation type implementation
NTL_vector_impl(Relation, vec_Relation)

NTL_io_vector_impl(Relation, vec_Relation)

NTL_eq_vector_impl(Relation, vec_Relation)


//vec_vec_Relation type declaration
NTL_vector_impl(vec_Relation, vec_vec_Relation)

NTL_eq_vector_impl(vec_Relation, vec_vec_Relation)

NTL_io_vector_impl(vec_Relation, vec_vec_Relation)

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


#include <NTL/new.h>

#include <pair_uintType_uintType.h>

//NTL_START_IMPL

NTL_pair_impl(uint_type, uint_type, pair_uintType_uintType)
NTL_pair_io_impl(uint_type, uint_type, pair_uintType_uintType)
NTL_pair_eq_impl(uint_type, uint_type, pair_uintType_uintType)

NTL_vector_impl(pair_uintType_uintType, vec_pair_uintType_uintType)
NTL_io_vector_impl(pair_uintType_uintType, vec_pair_uintType_uintType)
NTL_eq_vector_impl(pair_uintType_uintType, vec_pair_uintType_uintType)

//NTL_END_IMPL

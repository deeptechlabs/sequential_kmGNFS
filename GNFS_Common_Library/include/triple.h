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


#ifndef _NTL_LIKE_TRIPLE_H
#define	_NTL_LIKE_TRIPLE_H

#include <NTL/tools.h>

//Template like makros for triples, following NTL's method of pairs

#define NTL_triple_decl(S,T,U,triple_S_T_U)  \
class triple_S_T_U {  \
public:  \
   S a;  \
   T b;  \
   U c;  \
  \
   triple_S_T_U() { }  \
   triple_S_T_U(const triple_S_T_U& l__x) : a(l__x.a), b(l__x.b), c(l__x.c)  { } \
   triple_S_T_U& operator=(const triple_S_T_U& l__x) { a = l__x.a; b = l__x.b; c = l__x.c; return *this; } \
   triple_S_T_U(const S& l__x, const T& l__y, const U& l__z) : a(l__x), b(l__y), c(l__z) { }  \
   ~triple_S_T_U() { }  \
};  \
  \
inline triple_S_T_U cons(const S& l__x, const T& l__y, const U& l__z) { return triple_S_T_U(l__x, l__y, l__z); } \




#define NTL_triple_io_decl(S,T,U,triple_S_T_U) \
NTL_SNS istream& operator>>(NTL_SNS istream&, triple_S_T_U&);  \
  \
NTL_SNS ostream& operator<<(NTL_SNS ostream&, const triple_S_T_U&);  \



#define NTL_triple_eq_decl(S,T,U,triple_S_T_U)  \
inline long operator==(const triple_S_T_U& l__x, const triple_S_T_U& l__y)  \
   { return l__x.a == l__y.a && l__x.b == l__y.b && l__x.c == l__y.c; }  \
inline long operator!=(const triple_S_T_U& l__x, const triple_S_T_U& l__y) \
   { return !(l__x == l__y); }  \



// For compatability...
#define NTL_triple_impl(S,T,U,triple_S_T_U)


#define NTL_triple_io_impl(S,T,U,triple_S_T_U)  \
NTL_SNS istream& operator>>(NTL_SNS istream& l__s, triple_S_T_U& l__x)  \
{  \
   long l__c;  \
  \
   if (!l__s) NTL_NNS Error("bad pair input");  \
  \
   l__c = l__s.peek();  \
   while (NTL_NNS IsWhiteSpace(l__c)) {  \
      l__s.get();  \
      l__c = l__s.peek();  \
   }  \
  \
   if (l__c != '[')  \
      NTL_NNS Error("bad pair input");  \
  \
   l__s.get();  \
  \
   if (!(l__s >> l__x.a))   \
      NTL_NNS Error("bad pair input");  \
   if (!(l__s >> l__x.b))  \
      NTL_NNS Error("bad pair input");  \
   if (!(l__s >> l__x.c))  \
      NTL_NNS Error("bad pair input");  \
  \
   l__c = l__s.peek();  \
   while (NTL_NNS IsWhiteSpace(l__c)) {  \
      l__s.get();  \
      l__c = l__s.peek();  \
   }  \
  \
   if (l__c != ']')  \
      NTL_NNS Error("bad pair input");  \
  \
   l__s.get();  \
  \
   return l__s;  \
}  \
  \
NTL_SNS ostream& operator<<(NTL_SNS ostream& l__s, const triple_S_T_U& l__x)  \
{  \
   return l__s << "[" << l__x.a << " " << l__x.b << " " << l__x.c << "]";  \
}  \



// For compatability...
#define NTL_triple_eq_impl(S,T,U,triple_S_T_U)

#endif	/* _NTL_LIKE_TRIPLE_H */


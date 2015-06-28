#ifndef _MATH_ROUTINES_HPP
#define _MATH_ROUTINES_HPP

namespace nissa
{
  double lfact(double n);
  double metro_tresh(double arg);
  int metro_test(double arg);
  int factorize(int *list,int N);
  int log2N(int N);
  void matrix_determinant(complex d,complex *m,int n);
  int bitrev(int in,int l2n);
  int find_max_pow2(int a);

  template <class T> T sqr(T a){return a*a;}
  template <class T> T cube(T a){return a*a*a;};
}

#endif

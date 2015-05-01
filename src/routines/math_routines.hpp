#ifndef _NISSA_MATH_HPP
#define _NISSA_MATH_HPP

namespace nissa
{
  double lfact(double n);
  double metro_tresh(double arg);
  int metro_test(double arg);
  int factorize(int *list,int N);
  int log2N(int N);
  void matrix_determinant(complex d,complex *m,int n);

  template <class T> T sqr(T a){return a*a;}
  template <class T> T cube(T a){return a*a*a;};
}

#endif

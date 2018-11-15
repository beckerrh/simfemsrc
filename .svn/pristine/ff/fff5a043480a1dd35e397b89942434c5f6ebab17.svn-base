#ifndef __alat_alatmath_h
#define __alat_alatmath_h

#include  <string>
#include  <cmath>

/*--------------------------------------------------------------------------*/
namespace alatmath
{
  inline int binomial(int n, int k)
  {
    int res = 1;
    k = min(k, n-k);
    // if ( k > n - k ) {k = n - k;}
    for (int i = 0; i < k; ++i)
    {
      res *= (n - i);
      res /= (i + 1);
    }
    return res;
  }
  inline double phi(double x, double y)
  {
    if(y >= 0.0)
    {
      return atan2(y, x);
    }
    else
    {
      return 2.0*M_PI-atan2(-y, x);
    }
  }
}

/*--------------------------------------------------------------------------*/

#endif

#include  "Alat/random.hpp"
#include  <stdlib.h>

namespace alat
{
  double rangedRandomNumber(double fMin, double fMax)
  {
    double random =   (double)rand();
    return fMin + random * ( fMax - fMin )/ (double) RAND_MAX;
  }
}
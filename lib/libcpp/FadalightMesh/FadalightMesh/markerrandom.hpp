#ifndef __FadalightMesh_MarkerRandom_h
#define __FadalightMesh_MarkerRandom_h

#include  "markerbase.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class MarkerRandom : public MarkerBase
  {
private:
    int _n;
    void _mark();

public:
    MarkerRandom(int n, double marking_parameter) : MarkerBase(marking_parameter), _n(n) {}
    ~MarkerRandom() {}
    std::string getClassName() const
    {
      return "MarkerRandom";
    }

    void write(std::string outfilename, arma::file_type datatype = arma::arma_binary);
  };
}

/*--------------------------------------------------------------------------*/

#endif

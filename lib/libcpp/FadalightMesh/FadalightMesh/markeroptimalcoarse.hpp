#ifndef __FadalightMesh_MarkerOptimalCoarse_h
#define __FadalightMesh_MarkerOptimalCoarse_h

#include  "markerbase.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class MarkerOptimalCoarse : public MarkerBase
  {
private:
    double _s, _p, _alpha;
    int _dimension;
    void _mark();
    double _f(int i, int n, double total, double sum) const;

public:
    MarkerOptimalCoarse(double marking_parameter, int dimension) : MarkerBase(marking_parameter), _dimension(dimension), _s(marking_parameter) {}
    ~MarkerOptimalCoarse() {}
    std::string getClassName() const
    {
      return "MarkerOptimalCoarse";
    }

    void write(std::string outfilename, arma::file_type datatype = arma::arma_binary);
  };
}

/*--------------------------------------------------------------------------*/

#endif

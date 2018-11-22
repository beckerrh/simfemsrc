#ifndef __FadalightMesh_MarkerCoarsePercentage_h
#define __FadalightMesh_MarkerCoarsePercentage_h

#include  "markerbase.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class MarkerCoarsePercentage : public MarkerBase
  {
private:
    void _mark();

public:
    MarkerCoarsePercentage(double marking_parameter) : MarkerBase(marking_parameter) {}
    ~MarkerCoarsePercentage() {}
    std::string getClassName() const
    {
      return "MarkerCoarsePercentage";
    }
    void write(std::string outfilename, arma::file_type datatype = arma::arma_binary);
  };
}

/*--------------------------------------------------------------------------*/

#endif

#ifndef __FadalightMesh_MarkerPercentage_h
#define __FadalightMesh_MarkerPercentage_h

#include  "markerbase.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class MarkerPercentage : public MarkerBase
  {
private:
    void _mark();

public:
    MarkerPercentage(double marking_parameter) : MarkerBase(marking_parameter) {}
    ~MarkerPercentage() {}
    std::string getClassName() const
    {
      return "MarkerPercentage";
    }
    void write(std::string outfilename, arma::file_type datatype = arma::arma_binary);
  };
}

/*--------------------------------------------------------------------------*/

#endif

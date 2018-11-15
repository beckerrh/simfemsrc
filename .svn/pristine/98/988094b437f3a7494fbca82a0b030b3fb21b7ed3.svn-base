#ifndef __FadalightMesh_MarkerMean_h
#define __FadalightMesh_MarkerMean_h

#include  "markerbase.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class MarkerMean : public MarkerBase
  {
  private:
    void _mark();

  public:
    MarkerMean(double marking_parameter) : MarkerBase(marking_parameter) {}
    ~MarkerMean() {}
    std::string getClassName() const {return "MarkerMean";}
    void write(std::string outfilename, arma::file_type datatype = arma::arma_binary);
  };
}

/*--------------------------------------------------------------------------*/

#endif

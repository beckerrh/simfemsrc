#ifndef __FadalightMesh_MarkerBase_h
#define __FadalightMesh_MarkerBase_h

#include  "markerinterface.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class MarkerBase : public MarkerInterface
  {
protected:
    double _marking_parameter;
    alat::armavec _indicator;
    alat::armaivec _marked_cells;
    alat::armaivec _C;

protected:
    void _sort();

public:
    MarkerBase(double marking_parameter);
    ~MarkerBase() {}
    std::string getClassName() const
    {
      return "MarkerBase";
    }

    alat::armavec& getIndicator()
    {
      return _indicator;
    }

    void write(std::string outfilename, arma::file_type datatype = arma::arma_binary);
  };
}

/*--------------------------------------------------------------------------*/

#endif

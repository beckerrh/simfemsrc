#ifndef __FadalightMesh_MarkerInterface_h
#define __FadalightMesh_MarkerInterface_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class MarkerInterface : public alat::InterfaceBase
  {
protected:
    std::string getInterfaceName() const
    {
      return "MarkerInterface";
    }

public:
    MarkerInterface() : alat::InterfaceBase() {}
    ~MarkerInterface() {}
    virtual alat::armavec& getIndicator() = 0;
    virtual void write(std::string outfilename, arma::file_type datatype = arma::arma_binary)=0;
  };
}

/*--------------------------------------------------------------------------*/

#endif

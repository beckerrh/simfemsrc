#ifndef __FadalightMesh_CurvedBoundaryDescriptionInterface_h
#define __FadalightMesh_CurvedBoundaryDescriptionInterface_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/

namespace alat
{
  class Node;
}

namespace FadalightMesh
{
  class CurvedBoundaryDescriptionInterface : public alat::InterfaceBase
  {
protected:
    std::string getInterfaceName() const;

public:
    ~CurvedBoundaryDescriptionInterface();
    CurvedBoundaryDescriptionInterface();
    CurvedBoundaryDescriptionInterface( const CurvedBoundaryDescriptionInterface& curvedboundarydescriptioninterface);
    CurvedBoundaryDescriptionInterface& operator=( const CurvedBoundaryDescriptionInterface& curvedboundarydescriptioninterface);

    virtual double operator()(double x, double y, double z) const;
    virtual void grad(alat::Node& dst, const alat::Node& src) const = 0;
    virtual double newton(alat::Node&) const = 0;
    virtual double newton(alat::Node&, const alat::Node& x1, const alat::Node& x2) const=0;
    virtual void write(std::ostream& out, arma::file_type datatype = arma::arma_binary) const;
    virtual void read(std::istream& in);
    virtual alat::armavec& getParameters() = 0;
    virtual const alat::armavec& getParameters() const = 0;
    virtual CurvedBoundaryDescriptionInterface* clone() const = 0;
  };
}

/*--------------------------------------------------------------------------*/

#endif

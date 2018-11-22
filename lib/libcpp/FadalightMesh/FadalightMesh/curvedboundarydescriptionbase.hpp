#ifndef  __FadalightMesh_CurvedBoundaryDescriptionBase_h
#define  __FadalightMesh_CurvedBoundaryDescriptionBase_h

#include  "FadalightMesh/curvedboundarydescriptioninterface.hpp"

/*---------------------------------------------------*/

namespace FadalightMesh
{
  class CurvedBoundaryDescriptionBase : public FadalightMesh::CurvedBoundaryDescriptionInterface
  {
private:
    alat::armavec _parameters;

protected:
    double getParameter(int i) const;

public:
    ~CurvedBoundaryDescriptionBase();
    CurvedBoundaryDescriptionBase();
    CurvedBoundaryDescriptionBase(const CurvedBoundaryDescriptionBase& curvedboundarydescriptionbase);
    CurvedBoundaryDescriptionBase& operator=(const CurvedBoundaryDescriptionBase& curvedboundarydescriptionbase);

    alat::armavec& getParameters();
    const alat::armavec& getParameters() const;
    void grad(alat::Node& dst, const alat::Node& src) const;
    double newton(alat::Node&) const;
    double newton(alat::Node&, const alat::Node& x1, const alat::Node& x2) const;
    void write(std::ostream& out, arma::file_type datatype = arma::arma_binary) const;
    void read(std::istream& in);
  };
}

/*---------------------------------------------------*/

#endif

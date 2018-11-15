#ifndef  __FadalightMesh_CurvedBoundaryDescriptionPolynomial_h
#define  __FadalightMesh_CurvedBoundaryDescriptionPolynomial_h

#include  "Alat/node.hpp"
#include  "curvedboundarydescriptionbase.hpp"

/*---------------------------------------------------*/

namespace FadalightMesh
{
  class CurvedBoundaryDescriptionLinear2d : public CurvedBoundaryDescriptionBase
  {
    const alat::Node& _n1, _n2;
    double _a, _b, _c;
public:
    CurvedBoundaryDescriptionLinear2d(const alat::Node& n1,  const alat::Node& n2);
    CurvedBoundaryDescriptionLinear2d(const CurvedBoundaryDescriptionLinear2d& curvedboundarydescription);
    ~CurvedBoundaryDescriptionLinear2d();
    FadalightMesh::CurvedBoundaryDescriptionInterface* clone() const;
    std::string getClassName() const;
    double operator()(double x, double y, double z) const;
    void grad(alat::Node& dst, const alat::Node& src) const;
  };

  class CurvedBoundaryDescriptionQuadratic2d : public CurvedBoundaryDescriptionBase
  {
public:
    CurvedBoundaryDescriptionQuadratic2d();
    CurvedBoundaryDescriptionQuadratic2d(const CurvedBoundaryDescriptionQuadratic2d& curvedboundarydescription);
    ~CurvedBoundaryDescriptionQuadratic2d();
    FadalightMesh::CurvedBoundaryDescriptionInterface* clone() const;
    std::string getClassName() const;
    double operator()(double x, double y, double z) const;
    void grad(alat::Node& dst, const alat::Node& src) const;
  };
}

/*---------------------------------------------------*/

#endif

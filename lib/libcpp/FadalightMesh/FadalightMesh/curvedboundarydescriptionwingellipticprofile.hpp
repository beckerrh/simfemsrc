#ifndef  __FadalightMesh_CurvedBoundaryDescriptionWingEllipticProfil_h
#define  __FadalightMesh_CurvedBoundaryDescriptionWingEllipticProfil_h

#include  "curvedboundarydescriptionbase.hpp"

/*---------------------------------------------------*/

namespace FadalightMesh
{
  class CurvedBoundaryDescriptionWingEllipticProfil : public CurvedBoundaryDescriptionBase
  {
public:
    CurvedBoundaryDescriptionWingEllipticProfil();
    CurvedBoundaryDescriptionWingEllipticProfil(const CurvedBoundaryDescriptionWingEllipticProfil& curvedboundarydescription);
    ~CurvedBoundaryDescriptionWingEllipticProfil();
    FadalightMesh::CurvedBoundaryDescriptionInterface* clone() const;
    std::string getClassName() const;
    double operator()(double x, double y, double z) const;
    void grad(alat::Node& dst, const alat::Node& src) const;
  };
}

/*---------------------------------------------------*/

#endif

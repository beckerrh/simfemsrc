#ifndef  __FadalightMesh_CurvedBoundaryDescriptionNaca_h
#define  __FadalightMesh_CurvedBoundaryDescriptionNaca_h

#include  "curvedboundarydescriptionbase.hpp"

/*---------------------------------------------------*/

namespace FadalightMesh
{
  class CurvedBoundaryDescriptionNaca : public CurvedBoundaryDescriptionBase
  {
public:
    CurvedBoundaryDescriptionNaca();
    CurvedBoundaryDescriptionNaca(const CurvedBoundaryDescriptionNaca& curvedboundarydescription);
    ~CurvedBoundaryDescriptionNaca();
    FadalightMesh::CurvedBoundaryDescriptionInterface* clone() const;
    std::string getClassName() const;
    double operator()(double x, double y, double z) const;
    void grad(alat::Node& dst, const alat::Node& src) const;
  };
}

/*---------------------------------------------------*/

#endif

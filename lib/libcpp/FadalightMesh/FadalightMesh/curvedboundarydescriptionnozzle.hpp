#ifndef  __FadalightMesh_CurvedBoundaryDescriptionNozzle_h
#define  __FadalightMesh_CurvedBoundaryDescriptionNozzle_h

#include  "curvedboundarydescriptionbase.hpp"

/*---------------------------------------------------*/

namespace FadalightMesh
{
  class CurvedBoundaryDescriptionNozzle : public CurvedBoundaryDescriptionBase
  {
public:
    CurvedBoundaryDescriptionNozzle();
    CurvedBoundaryDescriptionNozzle(const CurvedBoundaryDescriptionNozzle& curvedboundarydescription);
    ~CurvedBoundaryDescriptionNozzle();
    FadalightMesh::CurvedBoundaryDescriptionInterface* clone() const;
    std::string getClassName() const;
    double operator()(double x, double y, double z) const;
    void grad(alat::Node& dst, const alat::Node& src) const;
  };
}

/*---------------------------------------------------*/

#endif

#include  "FadalightMesh/curvedboundarydescriptionconstructor.hpp"
#include  "FadalightMesh/curvedboundarydescriptionnaca.hpp"
#include  "FadalightMesh/curvedboundarydescriptionnozzle.hpp"
#include  "FadalightMesh/curvedboundarydescriptionpolynomial.hpp"
#include  "FadalightMesh/curvedboundarydescriptionwingellipticprofile.hpp"
#include  "Alat/stringvector.hpp"
#include  <cassert>

using namespace FadalightMesh;
using namespace std;

/*---------------------------------------------------*/
CurvedBoundaryDescriptionConstructor::~CurvedBoundaryDescriptionConstructor() {}
CurvedBoundaryDescriptionConstructor::CurvedBoundaryDescriptionConstructor() {}

FadalightMesh::CurvedBoundaryDescriptionInterface* CurvedBoundaryDescriptionConstructor::newDescription(std::string name) const
{
  if(name == "CurvedBoundaryDescriptionQuadratic2d")
  {
    return new CurvedBoundaryDescriptionQuadratic2d;
  }
  else if(name == "CurvedBoundaryDescriptionNozzle")
  {
    return new CurvedBoundaryDescriptionNozzle;
  }
  else if(name == "CurvedBoundaryDescriptionWingEllipticProfil")
  {
    return new CurvedBoundaryDescriptionWingEllipticProfil;
  }
  else if(name == "CurvedBoundaryDescriptionNaca")
  {
    return new CurvedBoundaryDescriptionNaca;
  }
  else
  {
    std::cerr << "*** CurvedBoundaryDescriptionConstructor() unkwnown CurvedBoundaryDescription " << name << "\n";
    assert(0);
  }
}
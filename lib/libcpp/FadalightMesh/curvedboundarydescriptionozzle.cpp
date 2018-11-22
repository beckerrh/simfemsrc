#include  "Alat/node.hpp"
#include  "FadalightMesh/curvedboundarydescriptionnozzle.hpp"
#include  <math.h>
#include  <cassert>

using namespace FadalightMesh;

/*---------------------------------------------------*/

CurvedBoundaryDescriptionNozzle::CurvedBoundaryDescriptionNozzle() : CurvedBoundaryDescriptionBase()
{}

CurvedBoundaryDescriptionNozzle::CurvedBoundaryDescriptionNozzle(const CurvedBoundaryDescriptionNozzle& curvedboundarydescription) : CurvedBoundaryDescriptionBase(curvedboundarydescription)
{
  assert(0);
}

CurvedBoundaryDescriptionNozzle::~CurvedBoundaryDescriptionNozzle() {}

FadalightMesh::CurvedBoundaryDescriptionInterface* CurvedBoundaryDescriptionNozzle::clone() const
{
  return new CurvedBoundaryDescriptionNozzle(*this);
}

std::string CurvedBoundaryDescriptionNozzle::getClassName() const
{
  return "CurvedBoundaryDescriptionNozzle";
}

/*---------------------------------------------------*/

double CurvedBoundaryDescriptionNozzle::operator()(double x, double y, double z) const
{
  const alat::armavec& p = getParameters();
  assert(p.size() == 1);

  // std::cerr << "CurvedBoundaryDescriptionNozzle::operator() x="<<x<<" y="<< y <<" p="<<p[0]<<" y-p="<<y-p[0]<<"\n";
  if( ( x >= 0.0 ) && ( x <= 4.0 ) )
  {
    return y - p[0]*0.25*( cos(0.5*M_PI*x)+3.0 );
  }
  return y-p[0];
}

/*---------------------------------------------------*/

void CurvedBoundaryDescriptionNozzle::grad(alat::Node& dst, const alat::Node& src) const
{
  const alat::armavec& p = getParameters();
  dst.y() = 1.0;
  dst.z() = 0.0;
  dst.x() = 0.0;
  double x = src.x();
  if( ( x >= 0.0 ) && ( x <= 4.0 ) )
  {
    dst.x() = -p[0]*0.125*M_PI*sin( 0.5*M_PI*x );
  }
}

#include  "FadalightMesh/curvedboundarydescriptioninterface.hpp"
#include  <cassert>

using namespace FadalightMesh;

/*--------------------------------------------------------------------------*/

CurvedBoundaryDescriptionInterface::~CurvedBoundaryDescriptionInterface()
{}

/*--------------------------------------------------------------------------*/

CurvedBoundaryDescriptionInterface::CurvedBoundaryDescriptionInterface() : alat::InterfaceBase()
{}

/*--------------------------------------------------------------------------*/

CurvedBoundaryDescriptionInterface::CurvedBoundaryDescriptionInterface( const CurvedBoundaryDescriptionInterface& curvedboundarydescriptioninterface) : alat::InterfaceBase(curvedboundarydescriptioninterface)
{}

/*--------------------------------------------------------------------------*/

CurvedBoundaryDescriptionInterface& CurvedBoundaryDescriptionInterface::operator=( const CurvedBoundaryDescriptionInterface& curvedboundarydescriptioninterface)
{
  InterfaceBase::operator=(curvedboundarydescriptioninterface);
  assert(0);
  return *this;
}

/*--------------------------------------------------------------------------*/

std::string CurvedBoundaryDescriptionInterface::getInterfaceName() const
{
  return "CurvedBoundaryDescriptionInterface";
}

/*--------------------------------------------------------------------------*/

double CurvedBoundaryDescriptionInterface::operator()(double x, double y, double z) const
{
  _notWritten("operator");
  return 0.0;
}

/*--------------------------------------------------------------------------*/


void CurvedBoundaryDescriptionInterface::write(std::ostream& out, arma::file_type datatype) const
{
  _notWritten("write");
}

/*--------------------------------------------------------------------------*/

void CurvedBoundaryDescriptionInterface::read(std::istream& in)
{
  _notWritten("read");
}

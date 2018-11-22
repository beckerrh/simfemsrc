#include  "Perulangan/vectorinterface.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
VectorInterface::~VectorInterface(){}
VectorInterface::VectorInterface() : alat::InterfaceBase(){}
VectorInterface::VectorInterface( const perulangan::VectorInterface& vectorinterface) : alat::InterfaceBase(vectorinterface){}
VectorInterface* VectorInterface::clone() const
{
  _notWritten("clone");
  return NULL;
}
std::string VectorInterface::getInterfaceName() const
{
  return "VectorInterface";
}
VectorInterface& VectorInterface::operator=( const perulangan::VectorInterface& vectorinterface)
{
  InterfaceBase::operator=(vectorinterface);
  assert(0);
  return *this;
}
/*--------------------------------------------------------------------------*/
double VectorInterface::norm() const{_notWritten("norm");return 0;}
double VectorInterface::scalarProduct(const perulangan::VectorInterface* v) const{_notWritten("scalarProduct");return 0;}
void VectorInterface::equal(const perulangan::VectorInterface* v){_notWritten("equal");}
void VectorInterface::equal(double d){_notWritten("equal");}
void VectorInterface::add(const double& d, const perulangan::VectorInterface* v){_notWritten("add");}
void VectorInterface::scale(const double& d){_notWritten("scale");}
void VectorInterface::fillzeros(){_notWritten("fillzeros");}

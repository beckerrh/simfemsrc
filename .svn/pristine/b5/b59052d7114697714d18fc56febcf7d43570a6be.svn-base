#include  "Solvers/integrationformulainterface.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
IntegrationFormulaInterface::~IntegrationFormulaInterface(){}
IntegrationFormulaInterface::IntegrationFormulaInterface(int n) : alat::InterfaceBase()
{
  _w.set_size(n);
  _c.set_size(n);
}
IntegrationFormulaInterface::IntegrationFormulaInterface( const IntegrationFormulaInterface& integrationformulainterface) : alat::InterfaceBase(integrationformulainterface)
{
  assert(0);
}
IntegrationFormulaInterface& IntegrationFormulaInterface::operator=( const IntegrationFormulaInterface& integrationformulainterface)
{
  InterfaceBase::operator=(integrationformulainterface);
  assert(0);
  return *this;
}
std::string IntegrationFormulaInterface::getInterfaceName() const
{
  return "solvers::IntegrationFormulaInterface";
}
std::string IntegrationFormulaInterface::getClassName() const
{
  return "IntegrationFormulaInterface";
}
/*--------------------------------------------------------------------------*/
int IntegrationFormulaInterface::n() const
{
  return _w.size();
}
double IntegrationFormulaInterface::weight(int k) const
{
  return _w[k];
}
const alat::Node& IntegrationFormulaInterface::point(int k)  const
{
  return _c[k];
}

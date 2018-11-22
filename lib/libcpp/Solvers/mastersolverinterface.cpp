#include  "Solvers/mastersolverinterface.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
MasterSolverInterface::~MasterSolverInterface() {}
MasterSolverInterface::MasterSolverInterface(): alat::InterfaceBase(){}
MasterSolverInterface::MasterSolverInterface( const MasterSolverInterface& mastersolverinterface): alat::InterfaceBase(mastersolverinterface)
{
assert(0);
}
MasterSolverInterface& MasterSolverInterface::operator=( const MasterSolverInterface& mastersolverinterface)
{
  assert(0);
  alat::InterfaceBase::operator=(mastersolverinterface);
  return *this;
}
std::string MasterSolverInterface::getClassName() const 
{
  return "MasterSolverInterface";
}
/*--------------------------------------------------------------------------*/

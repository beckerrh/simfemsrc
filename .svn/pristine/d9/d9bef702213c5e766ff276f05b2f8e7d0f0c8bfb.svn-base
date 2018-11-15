#include  "Perulangan/nonlinearsolvervisitorinterface.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
NonlinearSolverVisitorInterface::~NonlinearSolverVisitorInterface(){}
NonlinearSolverVisitorInterface::NonlinearSolverVisitorInterface() : alat::InterfaceBase(){}
NonlinearSolverVisitorInterface::NonlinearSolverVisitorInterface( const NonlinearSolverVisitorInterface& nonlinearsolvervisitorinterface) : alat::InterfaceBase(nonlinearsolvervisitorinterface)
{
  ( *this ).operator=(nonlinearsolvervisitorinterface);
}
NonlinearSolverVisitorInterface& NonlinearSolverVisitorInterface::operator=( const NonlinearSolverVisitorInterface& nonlinearsolvervisitorinterface)
{
  InterfaceBase::operator=(nonlinearsolvervisitorinterface);
  assert(0);
  return *this;
}
std::string NonlinearSolverVisitorInterface::getClassName() const
{
  return "NonlinearSolverVisitorInterface";
}
std::string NonlinearSolverVisitorInterface::getInterfaceName() const
{
  return "NonlinearSolverVisitorInterface";
}
NonlinearSolverVisitorInterface* NonlinearSolverVisitorInterface::clone() const
{
  assert(0);
  // return new NonlinearSolverVisitorInterface(*this);
}

/*--------------------------------------------------------------------------*/
std::ostream& NonlinearSolverVisitorInterface::printLoopInformation(std::ostream& os) const{_notWritten("printLoopInformation");return os;}
// std::string NonlinearSolverVisitorInterface::getVectorType() const{_notWritten("getVectorType");return "unknowns";}
void NonlinearSolverVisitorInterface::setLavrentievParameter(double parameter) const {_notWritten("setLavrentievParameter");}
double NonlinearSolverVisitorInterface::computeNormSquaredLavrientiev(perulanganEnums::residualstatus& status, const alat::GhostVector& u, const alat::GhostVector& du) const{ _notWritten("computeNormSquaredLavrientiev"); return 0.0;}
void NonlinearSolverVisitorInterface::computeLinearization(perulanganEnums::residualstatus& status, alat::GhostVector& h, const alat::GhostVector& u, const alat::GhostVector& du) const {_notWritten("computeLinearization");}

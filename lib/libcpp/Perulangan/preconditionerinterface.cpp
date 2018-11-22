#include  "Perulangan/preconditionerinterface.hpp"
#include  "Alat/stringvector.hpp"
#include  <cassert>
#include  <iostream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
PreconditionerInterface::~PreconditionerInterface(){}
PreconditionerInterface::PreconditionerInterface() : alat::InterfaceBase(){}
PreconditionerInterface::PreconditionerInterface( const PreconditionerInterface& vectorinterface) : alat::InterfaceBase(vectorinterface)
{
  assert(0);
}
PreconditionerInterface& PreconditionerInterface::operator=( const PreconditionerInterface& vectorinterface)
{
  InterfaceBase::operator=(vectorinterface);
  assert(0);
  return *this;
}
std::string PreconditionerInterface::getInterfaceName() const
{
  return "PreconditionerInterface";
}
std::string PreconditionerInterface::getClassName() const
{
  return "PreconditionerInterface";
}

std::ostream& PreconditionerInterface::printLoopInformation(std::ostream& os) const
{
  os << getClassName();
  return os;
}

perulangan::IterativeSolverVisitorInterface* PreconditionerInterface::getVisitor() {_notWritten("getVisitor");}
const perulangan::IterativeSolverVisitorInterface* PreconditionerInterface::getVisitor() const {_notWritten("getVisitor");}

/*--------------------------------------------------------------------------*/
// void PreconditionerInterface::basicInit(const alat::ParameterFile* parameterfile, std::string blockname, perulangan::IterativeSolverVisitorInterface* visitor)
// {
//   _notWritten("basicInit");
// }
//
void PreconditionerInterface::reInit()
{
  _notWritten("reInit");
}

void PreconditionerInterface::computePreconditioner()
{
  _notWritten("computePreconditioner");
}

void PreconditionerInterface::solveApproximate(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f, int iteration) const
{
  _notWritten("solveApproximate");
}

std::ostream& PreconditionerInterface::write(std::ostream& os) const
{
  _notWritten("write");
}

void PreconditionerInterface::fillzeros() const
{
  _notWritten("fillzeros");
}

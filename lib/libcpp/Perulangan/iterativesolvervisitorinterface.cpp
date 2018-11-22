#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  <cassert>
#include  <iostream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
IterativeSolverVisitorInterface::~IterativeSolverVisitorInterface(){}
IterativeSolverVisitorInterface::IterativeSolverVisitorInterface() : alat::InterfaceBase(){}
IterativeSolverVisitorInterface::IterativeSolverVisitorInterface( const IterativeSolverVisitorInterface& iterativesolvervisitorinterface) : alat::InterfaceBase(iterativesolvervisitorinterface)
{
  assert(0);
}

IterativeSolverVisitorInterface& IterativeSolverVisitorInterface::operator=( const IterativeSolverVisitorInterface& iterativesolvervisitorinterface)
{
  InterfaceBase::operator=(iterativesolvervisitorinterface);
  assert(0);
  return *this;
}

std::string IterativeSolverVisitorInterface::getInterfaceName() const
{
  return "IterativeSolverVisitorInterface";
}

std::string IterativeSolverVisitorInterface::getClassName() const
{
  return "IterativeSolverVisitorInterface";
}

IterativeSolverVisitorInterface* IterativeSolverVisitorInterface::clone() const
{
  assert(0);
  // return new IterativeSolverVisitorInterface(*this);
}

std::string IterativeSolverVisitorInterface::getVectorType() const
{
  return "unknowns";
}

// void IterativeSolverVisitorInterface::basicInit(const alat::ParameterFile* parameterfile, std::string blockname)
// {
//   _notWritten("basicInit");
// }

void IterativeSolverVisitorInterface::vectorEqual(alat::GhostVector& r, const alat::GhostVector& f) const
{
  _notWritten("vectorEqual");
}

void IterativeSolverVisitorInterface::vectorZero(alat::GhostVector& v) const
{
  _notWritten("vectorZero");
}

void IterativeSolverVisitorInterface::vectorAdd(alat::GhostVector& p, double d, const alat::GhostVector& q) const
{
  _notWritten("vectorAdd");
}

void IterativeSolverVisitorInterface::vectorScale(alat::GhostVector& r, double d) const
{
  _notWritten("vectorScale");
}

void IterativeSolverVisitorInterface::matrixVectorProduct(const alat::GhostMatrix& A, alat::GhostVector& r, const alat::GhostVector& u, double d) const
{
  _notWritten("matrixVectorProduct");
}

void IterativeSolverVisitorInterface::postProcess(alat::GhostVector& u) const
{}

/*-------------------------------------------------------------*/

double IterativeSolverVisitorInterface::vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const
{
  _notWritten("scalarProduct");
  return 0.0;
}

/*-------------------------------------------------------------*/

double IterativeSolverVisitorInterface::vectorNorm(const alat::GhostVector& r) const
{
  return sqrt( vectorDot(r, r) );
}

/*-------------------------------------------------------------*/

void IterativeSolverVisitorInterface::residual(const alat::GhostMatrix& A, alat::GhostVector& r, const alat::GhostVector& u, const alat::GhostVector& f) const
{
  vectorEqual(r, f);
  matrixVectorProduct(A, r, u, -1.0);
}

/*-------------------------------------------------------------*/

std::ostream& IterativeSolverVisitorInterface::printLoopInformation(std::ostream& os) const
{
  os << "\"" <<getClassName() <<"\"";
  return os;
}

// /*-------------------------------------------------------------*/
//
// void IterativeSolverVisitorInterface::matrixVectorProduct(perulangan::SystemVectorInterface* r, const perulangan::SystemVectorInterface* u, double d) const
// {
//   _notWritten("matrixVectorProduct(perulangan::SystemVectorInterface*)");
// }
//
// /*-------------------------------------------------------------*/
// const alat::armaivec& IterativeSolverVisitorInterface::getDomainsPermutation(int iteration) const
// {
//   _notWritten("getDomainsPermutation");
// }
//
// void IterativeSolverVisitorInterface::solveOnDomain(int idomain, const alat::GhostLinearSolver& linearsolverdomain, const alat::GhostMatrix& ghostmatrix, alat::GhostVector& u, const alat::GhostVector& f) const
// {
//   _notWritten("solveOnDomain");
// }
//
// void IterativeSolverVisitorInterface::vectorEqualOnDomain(int idomain, alat::GhostVector& u, const alat::GhostVector& f) const
// {
//   _notWritten("vectorEqualOnDomain");
// }
//
// void IterativeSolverVisitorInterface::matrixVectorProductCoupling(int i, const alat::GhostMatrix& ghostmatrix, alat::GhostVector& u, const alat::GhostVector& f, double d) const
// {
//   _notWritten("matrixVectorProductCoupling");
// }
//
// void IterativeSolverVisitorInterface::smoothInterface(int idomain, alat::GhostVector& u) const
// {
//   _notWritten("smoothInterface");
// }
//
// void IterativeSolverVisitorInterface::smoothInterfaceOnLevel(int level, alat::GhostVector& u) const
// {
//   _notWritten("smoothInterface");
// }

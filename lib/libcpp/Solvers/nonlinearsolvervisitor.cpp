#include  "Alat/ghostvector.hpp"
#include  "Solvers/nonlinearsolvervisitor.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
NonlinearSolverVisitor::~NonlinearSolverVisitor() {}
NonlinearSolverVisitor::NonlinearSolverVisitor(SolverInterface* solver): perulangan::NonlinearSolverVisitorInterface(), _solver(solver){}
NonlinearSolverVisitor::NonlinearSolverVisitor( const NonlinearSolverVisitor& nonlinearsolvervisitor): perulangan::NonlinearSolverVisitorInterface(nonlinearsolvervisitor)
{
assert(0);
}
NonlinearSolverVisitor& NonlinearSolverVisitor::operator=( const NonlinearSolverVisitor& nonlinearsolvervisitor)
{
  assert(0);
  perulangan::NonlinearSolverVisitorInterface::operator=(nonlinearsolvervisitor);
  return *this;
}
std::string NonlinearSolverVisitor::getClassName() const
{
  return "NonlinearSolverVisitor";
}
/*--------------------------------------------------------------------------*/
void NonlinearSolverVisitor::init()
{
}
/*--------------------------------------------------------------------------*/
std::string NonlinearSolverVisitor::getVectorType() const {return "solver";}
void NonlinearSolverVisitor::newVector(alat::GhostVector* u)
{
  _solver->enrolVector(alat::GhostVector(*u));
}
/*--------------------------------------------------------------------------*/
void NonlinearSolverVisitor::residual(perulanganEnums::residualstatus& status, alat::GhostVector& gr, const alat::GhostVector& gu, const alat::GhostVector& gf) const
{
  _solver->residual(status, gr, gu, gf);
}
int NonlinearSolverVisitor::solveLinear(perulanganEnums::iterationstatus& status, const alat::GhostLinearSolver& B, const alat::GhostMatrix& A, alat::GhostVector& du, const alat::GhostVector& r) const
{
  return _solver->solveLinear(status, B, A, du, r);
}
void NonlinearSolverVisitor::constructMatrixAndLinearSolvers(perulanganEnums::matrixstatus& status, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u)
{
  // return;
  _solver->constructMatrixAndLinearSolvers(status, B, A, u);
}
void NonlinearSolverVisitor::setLinearTolerance(double rtol, double gtol, alat::GhostLinearSolver& B)
{
  _solver->setLinearTolerance(rtol, gtol, B);
}
void NonlinearSolverVisitor::vectorZero(alat::GhostVector& gu) const
{
  _solver->vectorZero(gu);
}
void NonlinearSolverVisitor::vectorAdd(alat::GhostVector& gu, double s, const alat::GhostVector& gv) const
{
  _solver->vectorAdd(gu, s, gv);
}
void NonlinearSolverVisitor::vectorEqual(alat::GhostVector& gu, const alat::GhostVector& gv) const
{
  _solver->vectorEqual(gu, gv);
}
double NonlinearSolverVisitor::vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const
{
  return _solver->vectorDot(gu, gv);
}
std::ostream& NonlinearSolverVisitor::vectorWrite(std::ostream& os, const alat::GhostVector& r) const{assert(0); return os;}

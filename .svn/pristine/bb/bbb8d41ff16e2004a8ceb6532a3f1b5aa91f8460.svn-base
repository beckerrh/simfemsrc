#include  "Perulangan/simpleiterativesolver.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
SimpleIterativeSolver::~SimpleIterativeSolver()
{}

SimpleIterativeSolver::SimpleIterativeSolver() : IterativeSolverWithPreconditioner()
{}

SimpleIterativeSolver::SimpleIterativeSolver( const SimpleIterativeSolver& simpleiterativesolver) : IterativeSolverWithPreconditioner(simpleiterativesolver)
{
  assert(0);
}

SimpleIterativeSolver& SimpleIterativeSolver::operator=( const SimpleIterativeSolver& simpleiterativesolver)
{
  IterativeSolverWithPreconditioner::operator=(simpleiterativesolver);
  assert(0);
  return *this;
}

std::string SimpleIterativeSolver::getClassName() const
{
  return "SimpleIterativeSolver";
}

/*--------------------------------------------------------------------------*/

void SimpleIterativeSolver::solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const
{
  // std::cerr << "SimpleIterativeSolver::solve() u="<< u << " f="<< f<<  " prec="<<getPreconditioner()->getClassName()<<"\n";
  // std::cerr << "SimpleIterativeSolver::solve() precvisitor "<<getPreconditioner()->getVisitor()->getClassName()<<"\n";
  int iteration=0;
  getPreconditioner()->solveApproximate(status, A, u, f, iteration);
}

#include  "Solvers/mastersolver.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
MasterSolver::~MasterSolver() {}
MasterSolver::MasterSolver(): MasterSolverInterface(){}
MasterSolver::MasterSolver( const MasterSolver& mastersolver): MasterSolverInterface(mastersolver)
{
assert(0);
}
MasterSolver& MasterSolver::operator=( const MasterSolver& mastersolver)
{
  assert(0);
  MasterSolverInterface::operator=(mastersolver);
  return *this;
}
std::string MasterSolver::getClassName() const 
{
  return "MasterSolver";
}
/*--------------------------------------------------------------------------*/

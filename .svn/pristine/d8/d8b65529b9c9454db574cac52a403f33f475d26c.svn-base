#include  "Alat/ghostlinearsolver.hpp"
#include  <iostream>
#include  <cassert>

using namespace alat;

/*--------------------------------------------------------------------------*/
GhostLinearSolver::~GhostLinearSolver(){}
GhostLinearSolver::GhostLinearSolver(const std::string name, const std::string type, const alat::GhostMatrix matrix) : Ghost(name,type), _matrix(matrix) {}
GhostLinearSolver::GhostLinearSolver( const GhostLinearSolver& ghost) : Ghost(ghost), _matrix(ghost._matrix){}
GhostLinearSolver& GhostLinearSolver::operator=( const GhostLinearSolver& ghost)
{
  Ghost::operator=(ghost);
  assert(0);
  // _matrix = ghost._matrix;
  return *this;
}

std::string GhostLinearSolver::getClassName() const
{
  return "GhostLinearSolver";
}

std::ostream& alat::operator<<(std::ostream& os, const GhostLinearSolver& g)
{
  os << "(Name/type/matrix:) " << g.getName() <<"/"<< g.getType()<<"/"<< g.getMatrix();
  return os;
}
const alat::GhostMatrix& GhostLinearSolver::getMatrix() const {return _matrix;}

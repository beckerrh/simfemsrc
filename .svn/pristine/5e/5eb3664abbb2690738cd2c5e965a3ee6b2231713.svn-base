#include  "Alat/ghostvector.hpp"
#include  "Alat/vectorallvariables.hpp"
#include  "Perulangan/directlinearsolver.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
DirectLinearSolver::~DirectLinearSolver() {}
DirectLinearSolver::DirectLinearSolver(const alat::MatrixAllVariables* matrix, const solvers::MeshUnitWithDataInterface* solver): perulangan::LinearSolverInterface(), _matrixallvariables(matrix), _solver(solver){}
DirectLinearSolver::DirectLinearSolver( const DirectLinearSolver& directlinearsolver): perulangan::LinearSolverInterface(directlinearsolver)
{
assert(0);
}
DirectLinearSolver& DirectLinearSolver::operator=( const DirectLinearSolver& directlinearsolver)
{
  assert(0);
  perulangan::LinearSolverInterface::operator=(directlinearsolver);
  return *this;
}
std::string DirectLinearSolver::getClassName() const
{
  return "DirectLinearSolver";
}
/*--------------------------------------------------------------------------*/
bool DirectLinearSolver::hasIterationInfo() const {return false;}

/*--------------------------------------------------------------------------*/
void DirectLinearSolver::reInit()
{
  // std::cerr << "DirectLinearSolver::reInit()\n";
  _matrixallvariables->reInit(_matrixinone);
  _umfmatrix.reInit(&_matrixinone.matrix);
}
/*--------------------------------------------------------------------------*/
void DirectLinearSolver::compute()
{
  // std::cerr << "DirectLinearSolver::compute()\n";
  _matrixallvariables->compute(_matrixinone);
  _umfmatrix.computeLu();
}
/*--------------------------------------------------------------------------*/
void DirectLinearSolver::solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const
{
  alat::VectorAllVariables* out = _solver->getVector(u);
  const alat::VectorAllVariables* in = _solver->getVector(f);

  _matrixinone.in.fill(arma::fill::zeros);
  int n = out->size();
  assert(n==in->size());
  for(int i = 0; i < n; i++)
  {
    in->get(i)->addVectorRhsForDirectSolver(_matrixinone.offsets[i], _matrixinone.in);
  }
  // std::cerr << "MatrixAllVariables::solve() in=" << in<<"\n";
  // std::cerr << "MatrixAllVariables::solve() _in=" << _in.t();
  _umfmatrix.solve(_matrixinone.out, _matrixinone.in);
  // std::cerr << "MatrixAllVariables::solve() _out=" << _out.t();
  for(int i = 0; i < n; i++)
  {
    out->get(i)->setVectorFromDirectSolver(_matrixinone.offsets[i], _matrixinone.out);
  }
}

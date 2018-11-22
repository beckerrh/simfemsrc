#include  "Perulangan/newtonsimple.hpp"
#include  "Perulangan/nonlinearsolvervisitorinterface.hpp"
#include  <cassert>
#include  <cmath>
#include  <iomanip>
#include  <sstream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
NewtonSimple::~NewtonSimple(){}
NewtonSimple::NewtonSimple() : NonlinearSolver(){}
NewtonSimple::NewtonSimple( const NewtonSimple& newtonsimple) : NonlinearSolver(newtonsimple)
{
  assert(0);
}
NewtonSimple& NewtonSimple::operator=( const NewtonSimple& newtonsimple)
{
  NonlinearSolver::operator=(newtonsimple);
  assert(0);
  return *this;
}
std::string NewtonSimple::getClassName() const
{
  return "NewtonSimple";
}

NewtonSimple* NewtonSimple::clone() const
{
  return new NewtonSimple(*this);
}

/*--------------------------------------------------------------------------*/
int NewtonSimple::getNVectors() const
{
  return 2;
}

/*--------------------------------------------------------------------------*/
void NewtonSimple::solve(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f)
{
  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& du = getMemory(1);

  bool reached = 0;
  _info._printbufff = " dxnorm matrix (lintol nlin)\n";
  for(int iteration = 0; !reached; iteration++)
  {
    if( iteration == _info.getMaxiter() )
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusTooManyIterations;
      return;
    }
    getVisitor()->residual(newtonoutputdata.residual_status, r, u, f);
    _newtoninternaldata.residual =  sqrt(getVisitor()->vectorDot(r, r));
    if(newtonoutputdata.residual_status != perulanganEnums::ResidualStatusOk)
    {
      return;
    }
    _newtoninternaldata.iteration=iteration;
    if(_checkIteration( newtonoutputdata) )
    {
      return;
    }
    // std::cerr << "NewtonSimple::solve() _newtoninternaldata.newmatrix="<<_newtoninternaldata.newmatrix<<"\n";
    if(_newtoninternaldata.newmatrix)
    {
      getVisitor()->constructMatrixAndLinearSolvers(newtonoutputdata.matrix_status, B, A, u);
    }
    getVisitor()->setLinearTolerance(_newtoninternaldata.lintol, _info.getGlobalTol(), B);
    _newtoninternaldata.nliniter = getVisitor()->solveLinear(newtonoutputdata.linear_solver_status, B, A, du, r);
    newtonoutputdata.niter_linear += _newtoninternaldata.nliniter;
    _newtoninternaldata.dxnorm =  sqrt(getVisitor()->vectorDot(du, du));
    if( ( newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusMaxIter )||( newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusDiverged ) )
    {
      // std::cerr << "NewtonSimple()::solve() newtonoutputdata.linear_solver_status="<<perulanganEnums::iterationStatusToString(newtonoutputdata.linear_solver_status)<<"\n";
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusLinearNotOk;
      return;
    }
    getVisitor()->vectorAdd(u, 1.0, du);
  }
  assert(0);
}

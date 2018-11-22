#include  "Perulangan/newtonrb.hpp"
#include  "Perulangan/nonlinearsolvervisitorinterface.hpp"
#include  <cassert>
#include  <cmath>
#include  <iomanip>
#include  <sstream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
NewtonRb::~NewtonRb(){}
NewtonRb::NewtonRb(std::string linesearchtype) : NonlinearSolver(), _linesearchtype(linesearchtype)
{
  assert( ( linesearchtype == "monotonicty" )or(linesearchtype == "armijo") or(linesearchtype == "wolfe") );
}
NewtonRb::NewtonRb( const NewtonRb& newtonsimple) : NonlinearSolver(newtonsimple)
{
  assert(0);
}
NewtonRb& NewtonRb::operator=( const NewtonRb& newtonsimple)
{
  NonlinearSolver::operator=(newtonsimple);
  assert(0);
  return *this;
}
std::string NewtonRb::getClassName() const
{
  return "NewtonRb";
}
NewtonRb* NewtonRb::clone() const
{
  return new NewtonRb(*this);
}

/*--------------------------------------------------------------------------*/
int NewtonRb::getNVectors() const
{
  std::cerr << "NewtonRb::basicInit() _linesearchtype="<<_linesearchtype<<"\n";
  if(_linesearchtype == "monotonicty")
  {
    return 2;
  }
  else
  {
    return 3;
  }
}

/*--------------------------------------------------------------------------*/
void NewtonRb::solve(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f)
{
  _newtoninternaldata.use_linesearch=true;

  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& du = getMemory(1);

  bool reached = 0;
  double kappa;
  int irelax = 0;
  _info._printbufff = " dxnorm irelax matrix (lintol nlin)\n";
  // std::cerr << "FFFF NewtonRb::solve() \n";
  // getVisitor()->vectorWrite(std::cerr, f);
  // std::cerr << "\n";
  for(int iteration = 0; !reached; iteration++)
  {
    // std::cerr << "NewtonRb::solve() iteration ="<<iteration<< " u =" << getVisitor()->vectorNorm(u)<< " f =" << getVisitor()->vectorNorm(f) << "\n";

    if( iteration == _info.getMaxiter() )
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusTooManyIterations;
      return;
    }
    _newtoninternaldata.residual_old = _newtoninternaldata.residual;
    // std::cerr << "UUUU NewtonRb::solve() \n";
    // getVisitor()->vectorWrite(std::cerr, u);
    getVisitor()->residual(newtonoutputdata.residual_status, r, u, f);
    _newtoninternaldata.residual =  sqrt(getVisitor()->vectorDot(r, r));
    // std::cerr << "RRRR NewtonRb::solve() r= \n";
    // getVisitor()->vectorWrite(std::cerr, r);
    if(newtonoutputdata.residual_status != perulanganEnums::ResidualStatusOk)
    {
      return;
    }
    if(iteration)
    {
      if(_linesearchtype == "monotonicty")
      {
        _lineSearchMonotonicity(newtonoutputdata, u, du, r, f);
      }
      else if(_linesearchtype == "armijo")
      {
        _lineSearchArmijo(newtonoutputdata, u, du, r, f, kappa);
      }
      else
      {
        _lineSearchWolfe(newtonoutputdata, u, du, r, f, kappa);
      }
      // if(newtonoutputdata.newton_status == perulanganEnums::NewtonStatusMaxLineSearchAttained)
      // {
      //   return;
      // }
    }
    _newtoninternaldata.iteration=iteration;
    if( _checkIteration( newtonoutputdata) )
    {
      return;
    }

    // same ???
    if(_newtoninternaldata.newmatrix)
    {
      getVisitor()->constructMatrixAndLinearSolvers(newtonoutputdata.matrix_status, B, A, u);
    }
    getVisitor()->setLinearTolerance(_newtoninternaldata.lintol, _info.getGlobalTol(), B);
    _newtoninternaldata.nliniter = getVisitor()->solveLinear(newtonoutputdata.linear_solver_status, B, A, du, r);
    newtonoutputdata.niter_linear += _newtoninternaldata.nliniter;
    // _dxnorm =  getVisitor()->vectorNorm(du);
    _newtoninternaldata.dxnorm =  sqrt(getVisitor()->vectorDot(du, du));
    // std::cerr << "NewtonRb: LINEAR status=" << perulanganEnums::iterationStatusToString(newtonoutputdata.linear_solver_status) << "\n";
    if( ( newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusMaxIter )||( newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusDiverged )||( newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusProblem ) )
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusLinearNotOk;
      return;
    }
    // same ???



    if(_linesearchtype != "monotonicty")
    {
      alat::GhostVector& h = getMemory(2);
      perulanganEnums::residualstatus status_linearization = perulanganEnums::ResidualStatusOk;
      getVisitor()->vectorZero(h);
      getVisitor()->computeLinearization(status_linearization, h, u, du);
      if(status_linearization != perulanganEnums::ResidualStatusOk)
      {
        _error_string( "solve", "linearization not ok", perulanganEnums::residualStatusToString(status_linearization) );
      }
      kappa = getVisitor()->vectorDot(r, h);
    }
    getVisitor()->vectorAdd(u, 1.0, du);
  }
  assert(0);
}

#include  "Perulangan/newtonlinesearch.hpp"
#include  "Perulangan/nonlinearsolvervisitorinterface.hpp"
#include  <cassert>
#include  <cmath>
#include  <iomanip>
#include  <limits>
#include  <sstream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
NewtonLineSearch::~NewtonLineSearch(){}
NewtonLineSearch::NewtonLineSearch(std::string linesearchtype) : NonlinearSolver(), _linesearchtype(linesearchtype)
{
  assert( ( linesearchtype == "monotonicty" )or(linesearchtype == "armijo") or(linesearchtype == "wolfe") );
}
NewtonLineSearch::NewtonLineSearch( const NewtonLineSearch& newtonsimple) : NonlinearSolver(newtonsimple)
{
  assert(0);
}
NewtonLineSearch& NewtonLineSearch::operator=( const NewtonLineSearch& newtonsimple)
{
  NonlinearSolver::operator=(newtonsimple);
  assert(0);
  return *this;
}

std::string NewtonLineSearch::getClassName() const
{
  return "NewtonLineSearch";
}

NewtonLineSearch* NewtonLineSearch::clone() const
{
  return new NewtonLineSearch(*this);
}

/*--------------------------------------------------------------------------*/
int NewtonLineSearch::getNVectors() const
{
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
void NewtonLineSearch::solve(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f)
{
  _newtoninternaldata.use_linesearch = true;

  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& du = getMemory(1);

  bool reached = 0;
  double kappa;
  int irelax = 0;
  _info._printbufff = " dxnorm irelax matrix (lintol nlin)\n";
  // std::cerr << "FFFF NewtonLineSearch::solve() \n";
  // getVisitor()->vectorWrite(std::cerr, f);
  // std::cerr << "\n";
  std::cerr << "NewtonLineSearch::solve() info =" << _info << "\n";
  for(int iteration = 0; !reached; iteration++)
  {
    // std::cerr << "NewtonLineSearch::solve() iteration ="<<iteration<< " u =" << getVisitor()->vectorNorm(u)<< " f =" << getVisitor()->vectorNorm(f) << "\n";

    if( iteration == _info.getMaxiter() )
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusTooManyIterations;
      return;
    }
    _newtoninternaldata.residual_old = _newtoninternaldata.residual;
    // std::cerr << "UUUU NewtonLineSearch::solve() \n";
    // getVisitor()->vectorWrite(std::cerr, u);
    getVisitor()->residual(newtonoutputdata.residual_status, r, u, f);
    _newtoninternaldata.residual =  sqrt(getVisitor()->vectorDot(r, r));
    if( _newtoninternaldata.residual >= std::numeric_limits<double>::max() )
    {
      std::cerr << "#NewtonLineSearch::solve() _newtoninternaldata.residual=" << _newtoninternaldata.residual<<"\n";
      newtonoutputdata.residual_status = perulanganEnums::ResidualStatusNotOk;
    }
    // std::cerr << "RRRR NewtonLineSearch::solve() r= \n";
    // getVisitor()->vectorWrite(std::cerr, r);
    if(newtonoutputdata.residual_status != perulanganEnums::ResidualStatusOk)
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusResidualStatusNotOk;
      std::cerr << "NewtonLineSearch::solve() " << perulanganEnums::residualStatusToString( newtonoutputdata.residual_status )<<"\n";
      std::cerr << "NewtonLineSearch::solve() " << perulanganEnums::newtonStatusToString( newtonoutputdata.newton_status )<<"\n";
      if(iteration)
      {
        getVisitor()->vectorAdd(u, -1.0, du);
      }
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
      if(newtonoutputdata.residual_status != perulanganEnums::ResidualStatusOk)
      {
        newtonoutputdata.newton_status = perulanganEnums::NewtonStatusResidualStatusNotOk;
        std::cerr << "#NewtonLineSearch::solve() " << perulanganEnums::newtonStatusToString( newtonoutputdata.newton_status )<<"\n";
        assert(0);
        return;
      }
      // if(newtonoutputdata.newton_status == perulanganEnums::NewtonStatusMaxLineSearchAttained)
      // {
      //   return;
      // }
    }
    _newtoninternaldata.iteration = iteration;
    if( _checkIteration( newtonoutputdata) )
    {
      return;
    }

    // same ???
    if(_newtoninternaldata.newmatrix)
    {
      getVisitor()->constructMatrixAndLinearSolvers(newtonoutputdata.matrix_status, B, A, u);
    }
    if(newtonoutputdata.matrix_status == perulanganEnums::MatrixStatusNotOk)
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusResidualStatusNotOk;
      return;
    }
    getVisitor()->setLinearTolerance(_newtoninternaldata.lintol, _info.getGlobalTol(), B);
    // std::cerr << "NewtonLineSearch::solve() " << getVisitor()->getClassName() << "\n";
    _newtoninternaldata.nliniter = getVisitor()->solveLinear(newtonoutputdata.linear_solver_status, B, A, du, r);
    // std::cerr << "NewtonLineSearch::solve() _nliniter="<<_nliniter<<"\n";
    newtonoutputdata.niter_linear += _newtoninternaldata.nliniter;
    // _dxnorm =  getVisitor()->vectorNorm(du);
    _newtoninternaldata.dxnorm =  sqrt(getVisitor()->vectorDot(du, du));
    // std::cerr << "NewtonLineSearch: LINEAR status=" << perulanganEnums::iterationStatusToString(newtonoutputdata.linear_solver_status) << "\n";
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

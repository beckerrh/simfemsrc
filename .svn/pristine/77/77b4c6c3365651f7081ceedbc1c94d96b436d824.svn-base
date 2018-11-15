#include  "Perulangan/newtonlavrentievstab.hpp"
#include  "Perulangan/nonlinearsolvervisitorinterface.hpp"
#include  <cassert>
#include  <cmath>
#include  <iomanip>
#include  <sstream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
NewtonLavrentievStab::~NewtonLavrentievStab(){}
NewtonLavrentievStab::NewtonLavrentievStab() : NewtonLavrentiev(){}
NewtonLavrentievStab::NewtonLavrentievStab( const NewtonLavrentievStab& newtonsimple) : NewtonLavrentiev(newtonsimple)
{
  assert(0);
}
NewtonLavrentievStab& NewtonLavrentievStab::operator=( const NewtonLavrentievStab& newtonsimple)
{
  NewtonLavrentiev::operator=(newtonsimple);
  assert(0);
  return *this;
}
std::string NewtonLavrentievStab::getClassName() const
{
  return "NewtonLavrentievStab";
}
NewtonLavrentievStab* NewtonLavrentievStab::clone() const
{
  return new NewtonLavrentievStab(*this);
}

/*--------------------------------------------------------------------------*/
int NewtonLavrentievStab::getNVectors() const
{
  return 4;
}

/*--------------------------------------------------------------------------*/
bool NewtonLavrentievStab::_changeRegularization(int iteration, double& regularization, double modelfactor)
{
  bool parameterchanged = true;
  double rfactor = fabs(modelfactor-1.0);
  if( rfactor < 0.001)
  {
    regularization *= 0.1;
  }
  else if( rfactor < 0.1)
  {
    regularization *= 0.3;
  }
  else if( rfactor < 0.2)
  {
    regularization *= 0.4;
  }
  else if( rfactor < 0.3)
  {
    regularization *= 0.8;
  }
  else if( rfactor > 0.9)
  {
    regularization *= 20.0;
  }
  else if( rfactor > 0.8)
  {
    regularization *= 10.0;
  }
  else if( rfactor > 0.6)
  {
    regularization *= 5.0;
  }
  else if( rfactor > 0.5)
  {
    regularization *= 3.0;
  }
  else if( rfactor > 0.4)
  {
    regularization *= 2.0;
  }
  else if( rfactor > 0.3)
  {
    regularization *= 1.5;
  }
  else
  {
    parameterchanged = false;
  }
  return parameterchanged;
}

/*--------------------------------------------------------------------------*/
void NewtonLavrentievStab::solve(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f)
{
  _newtoninternaldata.use_linesearch = true;
  _newtoninternaldata.use_regularization = true;

  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& du = getMemory(1);
  alat::GhostVector& h = getMemory(2);
  alat::GhostVector& uold = getMemory(3);

  _newtoninternaldata.regularization = _firstregularization;
  bool reached = false, parameterchanged;
  double kappa, wnorm, respred2, rescor, rrpredicted;
  _info._printbufff = " dxnorm irelax [mf param] matrix (lintol nlin)\n";
  int irelax = 0;
  for(int iteration = 0; !reached; iteration++)
  {
    if( iteration == _info.getMaxiter() )
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusTooManyIterations;
      return;
    }
    _newtoninternaldata.residual_old = _newtoninternaldata.residual;
    getVisitor()->residual(newtonoutputdata.residual_status, r, u, f);
    _newtoninternaldata.residual =  sqrt(getVisitor()->vectorDot(r, r));
    if(newtonoutputdata.residual_status != perulanganEnums::ResidualStatusOk)
    {
      return;
    }
    getVisitor()->vectorZero(h);
    getVisitor()->vectorAdd(h, -1.0, r);

    // rescor = getVisitor()->vectorNorm(r);
    rescor =  sqrt(getVisitor()->vectorDot(r, r));

    if(iteration)
    {
      _lineSearchMonotonicity(newtonoutputdata, u, du, r, f);
    }
    if(newtonoutputdata.newton_status == perulanganEnums::NewtonStatusMaxLineSearchAttained)
    {
      getVisitor()->vectorEqual(u, uold);
      _error_string("solve", "maximum linesearch attained");
    }
    _newtoninternaldata.iteration = iteration;
    if( _checkIteration( newtonoutputdata) )
    {
      return;
    }
    if(iteration)
    {
      // rescor = getVisitor()->vectorNorm(r);
      // std::cerr << "_omega="<<_omega << "rescor="<<rescor<<" respred2="<<respred2<< " ?" << respred2/rescor<< "\n";
      rrpredicted = respred2;
      _newtoninternaldata.modelfactor = rrpredicted/rescor;

      // std::cerr << "_newtoninternaldata.modelfactor="<<_newtoninternaldata.modelfactor<<" _newtoninternaldata.regularization="<<_newtoninternaldata.regularization<<" rrpredicted="<<rrpredicted<<" rescor="<<rescor<<"\n";
      if(  newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusConverged )
      {
        parameterchanged = _changeRegularization(iteration, _newtoninternaldata.regularization, _newtoninternaldata.modelfactor);
      }
      if (_newtoninternaldata.regularization>_maxregularization)
      {
        std::cerr << "NewtonLavrentievStab::solve() maxregularization attained ("<<_maxregularization<<")\n";
        newtonoutputdata.newton_status = perulanganEnums::NewtonStatusDiverged;
        return;
      }
    }
    // std::cerr  <<"_newmatrix="<<_newmatrix<<" parameterchanged="<<parameterchanged<<" _newtoninternaldata.regularization"<<_newtoninternaldata.regularization<<"\n";
    if(_newtoninternaldata.newmatrix or parameterchanged)
    {
      getVisitor()->setLavrentievParameter(_newtoninternaldata.regularization);
      getVisitor()->constructMatrixAndLinearSolvers(newtonoutputdata.matrix_status, B, A, u);
    }
    // if(iteration)
    // {
    //   getVisitor()->integrateTimeRhs(r, u, du, 0.5);
    // }
    getVisitor()->setLinearTolerance(_newtoninternaldata.lintol, _info.getGlobalTol(), B);

    newtonoutputdata.linear_solver_status = perulanganEnums::IterationStatusRunning;
    while(newtonoutputdata.linear_solver_status != perulanganEnums::IterationStatusConverged)
    {
      _newtoninternaldata.nliniter = getVisitor()->solveLinear(newtonoutputdata.linear_solver_status, B, A, du, r);
      if(newtonoutputdata.linear_solver_status != perulanganEnums::IterationStatusConverged)
      {
        _newtoninternaldata.regularization *= 4.0;
        getVisitor()->setLavrentievParameter(_newtoninternaldata.regularization);
        getVisitor()->constructMatrixAndLinearSolvers(newtonoutputdata.matrix_status, B, A, u);
        std::cerr << "INCREASING DUE TO LINEAR SOLVER status"<< perulanganEnums::iterationStatusToString( newtonoutputdata.linear_solver_status ) << " regularization="<<_newtoninternaldata.regularization<<"\n";
      }
      if(newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusMaxIter)
      {
        break;
      }
    }


    newtonoutputdata.niter_linear += _newtoninternaldata.nliniter;
    // _dxnorm =  getVisitor()->vectorNorm(du);
    _newtoninternaldata.dxnorm =  sqrt(getVisitor()->vectorDot(du, du));
    // if( ( newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusMaxIter )||( newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusDiverged ) )
    if( newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusDiverged )
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusLinearNotOk;
      return;
    }
    perulanganEnums::residualstatus status;
    getVisitor()->computeLinearization(status, h, u, du);
    // respred2 =  getVisitor()->vectorNorm(h);
    respred2 =  sqrt(getVisitor()->vectorDot(h, h));

    getVisitor()->vectorEqual(uold, u);
    getVisitor()->vectorAdd(u, 1.0, du);
  }

  assert(0);
}

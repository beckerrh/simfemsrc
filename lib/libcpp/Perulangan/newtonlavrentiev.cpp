#include  "Perulangan/newtonlavrentiev.hpp"
#include  "Perulangan/nonlinearsolvervisitorinterface.hpp"
#include  <cassert>
#include  <cmath>
#include  <iomanip>
#include  <sstream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
NewtonLavrentiev::~NewtonLavrentiev(){}
NewtonLavrentiev::NewtonLavrentiev(std::string type) : NonlinearSolver(), _type(type)
{
  if(type == "test")
  {
    _LINEARIZATION_TEST = true;
    _USE_LINEARIZATION = true;
  }
  else if(type == "lin")
  {
    _LINEARIZATION_TEST = false;
    _USE_LINEARIZATION = true;
  }
  else if(type == "mat")
  {
    _LINEARIZATION_TEST = false;
    _USE_LINEARIZATION = false;
  }
  else
  {
    _error_string("NewtonLavrentiev", "unknown type \"" + type + "\"");
  }
}
NewtonLavrentiev::NewtonLavrentiev( const NewtonLavrentiev& newtonsimple) : NonlinearSolver(newtonsimple)
{
  assert(0);
}
NewtonLavrentiev& NewtonLavrentiev::operator=( const NewtonLavrentiev& newtonsimple)
{
  NonlinearSolver::operator=(newtonsimple);
  assert(0);
  return *this;
}
std::string NewtonLavrentiev::getClassName() const
{
  return "NewtonLavrentiev";
}
NewtonLavrentiev* NewtonLavrentiev::clone() const
{
  return new NewtonLavrentiev(*this);
}

/*--------------------------------------------------------------------------*/
int NewtonLavrentiev::getNVectors() const
{
  return 3;
}
//
// /*--------------------------------------------------------------------------*/
// void NewtonLavrentiev::basicInit(const alat::ParameterFile* parameterfile, std::string blockname)
// {
//   // std::cerr << "NewtonLavrentiev::basicInit() blockname="<<blockname<<"\n";
//   NonlinearSolver::basicInit(parameterfile, blockname);
//   alat::DataFormatHandler dataformathandler;
//   dataformathandler.insert("firstregularization", &_firstregularization, 1.0);
//   dataformathandler.insert("maxregularization", &_maxregularization, 1.e6);
//   alat::FileScanner filescanner(dataformathandler, parameterfile, blockname, 0);
// }

/*--------------------------------------------------------------------------*/
bool NewtonLavrentiev::_changeRegularization(int iteration, double& regularization, double modelfactor)
{
  double rfactor = fabs(modelfactor-1.0);
  // std::cerr << "_changeRegularization() _newtoninternaldata.regularization="<<_newtoninternaldata.regularization<<" modelfactor="<<modelfactor<<" rfactor="<<rfactor<<"\n";
  bool parameterchanged = true;
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
    regularization *= 100.0;
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
void NewtonLavrentiev::solve(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f)
{
  _newtoninternaldata.use_regularization=true;

  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& du = getMemory(1);
  alat::GhostVector& h = getMemory(2);

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
    if(_USE_LINEARIZATION)
    {
      getVisitor()->vectorZero(h);
      getVisitor()->vectorAdd(h, -1.0, r);
    }
    _newtoninternaldata.iteration=iteration;
    if(_checkIteration( newtonoutputdata) )
    {
      return;
    }
    if(iteration)
    {
      if(_LINEARIZATION_TEST or not _USE_LINEARIZATION)
      {
        rescor = fabs( getVisitor()->vectorDot(r, du) );
      }
      else
      {
        // rescor = getVisitor()->vectorNorm(r);
        rescor =  sqrt(getVisitor()->vectorDot(r, r));
      }
      if(not _USE_LINEARIZATION or _LINEARIZATION_TEST)
      {
        wnorm = getVisitor()->computeNormSquaredLavrientiev(newtonoutputdata.residual_status, u, du);
      }
      if(_LINEARIZATION_TEST)
      {
        // std::cerr << "respred2="<<respred2 << " =?= _newtoninternaldata.regularization*wnorm=" << _newtoninternaldata.regularization*wnorm << "\n";
        if( fabs(respred2-_newtoninternaldata.regularization*wnorm) > fmax(1e-4*respred2, 1e-8) )
        {
          std::cerr << "fabs(respred2-_newtoninternaldata.regularization*wnorm)="<<fabs(respred2-_newtoninternaldata.regularization*wnorm) << "\n";
          assert(0);
        }
      }
      if(_USE_LINEARIZATION)
      {
        // std::cerr << "_omega="<<_omega << "rescor="<<rescor<<" respred2="<<respred2<< " ?" << respred2/rescor<< "\n";
        rrpredicted = respred2;
      }
      else
      {
        rrpredicted = _newtoninternaldata.regularization*wnorm;
      }
      _newtoninternaldata.modelfactor = rrpredicted/rescor;

      // std::cerr << "_newtoninternaldata.modelfactor="<<_newtoninternaldata.modelfactor<<" _newtoninternaldata.regularization="<<_newtoninternaldata.regularization<<" rrpredicted="<<rrpredicted<<" rescor="<<rescor<<"\n";
      parameterchanged = _changeRegularization(iteration, _newtoninternaldata.regularization, _newtoninternaldata.modelfactor);

      if(_newtoninternaldata.regularization > 1e10)
      {
        std::cerr << "_newtoninternaldata.regularization wnorm " << _newtoninternaldata.regularization << " " << wnorm << "\n";
        std::cerr << "_newtoninternaldata.modelfactor rrpredicted " << _newtoninternaldata.modelfactor << " " << rrpredicted << "\n";
        assert(0);
      }
    }
    // std::cerr  <<"_newmatrix="<<_newmatrix<<" parameterchanged="<<parameterchanged<<" _newtoninternaldata.regularization"<<_newtoninternaldata.regularization<<"\n";
    if(_newtoninternaldata.newmatrix or parameterchanged)
    {
      getVisitor()->setLavrentievParameter(_newtoninternaldata.regularization);
      getVisitor()->constructMatrixAndLinearSolvers(newtonoutputdata.matrix_status, B, A, u);
    }
    getVisitor()->setLinearTolerance(_newtoninternaldata.lintol, _info.getGlobalTol(), B);
    _newtoninternaldata.nliniter = getVisitor()->solveLinear(newtonoutputdata.linear_solver_status, B, A, du, r);
    newtonoutputdata.niter_linear += _newtoninternaldata.nliniter;
    // _dxnorm =  getVisitor()->vectorNorm(du);
    _newtoninternaldata.dxnorm =  sqrt(getVisitor()->vectorDot(du, du));
    if( ( newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusMaxIter )||( newtonoutputdata.linear_solver_status == perulanganEnums::IterationStatusDiverged ) )
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusLinearNotOk;
      return;
    }
    if(_USE_LINEARIZATION)
    {
      perulanganEnums::residualstatus status;
      getVisitor()->computeLinearization(status, h, u, du);
      if(_LINEARIZATION_TEST)
      {
        respred2 =  fabs( getVisitor()->vectorDot(h, du) );
      }
      else
      {
        // respred2 =  getVisitor()->vectorNorm(h);
        respred2 =  sqrt(getVisitor()->vectorDot(h, h));
      }
    }
    getVisitor()->vectorAdd(u, 1.0, du);
  }

  assert(0);
}

#include  "Perulangan/linearsolverinterface.hpp"
#include  "Perulangan/nonlinearsolver.hpp"
#include  "Perulangan/nonlinearsolvervisitorinterface.hpp"
#include  <cassert>
#include  <cmath>
#include  <iomanip>
#include  <sstream>
#include  <limits>

using namespace perulangan;

int NonlinearSolver::_totalids = 0;

/*--------------------------------------------------------------------------*/
NewtonInternalData::NewtonInternalData()
{
  perulanganEnums::iterationstatus status_iteration = perulanganEnums::IterationStatusNone;
  newmatrix = false;
  use_linesearch = false;
  use_regularization = false;
  nliniter =-1; iteration = 0;
  irelax = 0;
  lintol = 0.0;
  dxnorm = 0.0;
  omega = 0.0;
  residual = 0.0;
  residual_old = 0.0;
  modelfactor = 0.0;
  regularization = 0.0;
}
std::ostream& perulangan::operator<<(std::ostream& os, const NewtonInternalData& nd)
{
  os << "newmatrix/lintol/dxnorm " <<nd.newmatrix<<"/"<<nd.lintol<<"/"<<nd.dxnorm;
  return os;
}


/*--------------------------------------------------------------------------*/
NonlinearSolver::~NonlinearSolver(){}
NonlinearSolver::NonlinearSolver() : NonlinearSolverInterface(), _visitor(NULL), _basicinitcalled(false), _info()
{
  _id = _totalids;
  _totalids++;
}
NonlinearSolver::NonlinearSolver( const NonlinearSolver& nonlinearsolver) : NonlinearSolverInterface(nonlinearsolver)
{
  ( *this ).operator=(nonlinearsolver);
}
NonlinearSolver& NonlinearSolver::operator=( const NonlinearSolver& nonlinearsolver)
{
  NonlinearSolverInterface::operator=(nonlinearsolver);
  assert(0);
  return *this;
}
std::string NonlinearSolver::getClassName() const
{
  return "NonlinearSolver";
}
NonlinearSolver* NonlinearSolver::clone() const
{
  assert(0);
  // return new NonlinearSolver(*this);
}
/*--------------------------------------------------------------------------*/
void NonlinearSolver::_setNewtonInputData()
{
  _info.getIterationInfoData().rtol = _newtoninputdata.rtol;
  _info.getIterationInfoData().gtol = _newtoninputdata.gtol;
  int maxiter = log(_newtoninputdata.rtol)/log(_newtoninputdata.rhomatrix);
  if(maxiter>_newtoninputdata.maxiter)
  {
    std::cerr << "***Changing maxiter " << _newtoninputdata.maxiter << "  -->  " << maxiter<<"\n";
    assert(0);
    _newtoninputdata.rhomatrix = maxiter;
  }
  _info.getIterationInfoData().maxiter = _newtoninputdata.maxiter;
  _info.getIterationInfoData().printstep = _newtoninputdata.printstep;
}
void NonlinearSolver::setNewtonInputData(const NewtonInputData newtoninputdata)
{
  _newtoninputdata=newtoninputdata;
  _setNewtonInputData();
}
const NewtonInputData NonlinearSolver::getNewtonInputData() const {return _newtoninputdata;}

/*--------------------------------------------------------------------------*/
void NonlinearSolver::init()
{
  assert(_basicinitcalled == false);
  // visiteur
  assert(_visitor);
  getVisitor()->init();
  // info
  // getIterationInfo().basicInit(parameterfile, blockname);
  _setNewtonInputData();
	_info.setId(getClassName());
  memory();
  _basicinitcalled = true;
  // std::cerr << "NonlinearSolver::basicInit() _monotonyfactor="<<_monotonyfactor<<" _maxnlinesearch="<<_maxnlinesearch<<"\n";
  // assert(0);
}

/*--------------------------------------------------------------------------*/
std::ostream& NonlinearSolver::printLoopInformation(std::ostream& os) const
{
  os << "\"" << getClassName() << "\" ";
  getVisitor()->printLoopInformation(os);
  os << " info: ";
  _info.printLoopInformation(os);
  os << " ";
  return os;
}

/*--------------------------------------------------------------------------*/
void NonlinearSolver::setVisitorPointer(std::shared_ptr<perulangan::NonlinearSolverVisitorInterface> visitor) {_visitor=visitor;}

const perulangan::NonlinearSolverVisitorInterface* NonlinearSolver::getVisitor() const
{
  assert(_visitor);
  return _visitor.get();
}
perulangan::NonlinearSolverVisitorInterface* NonlinearSolver::getVisitor()
{
  assert(_visitor);
  return _visitor.get();
}

/*--------------------------------------------------------------------------*/
alat::GhostVector& NonlinearSolver::getMemory(int i) const
{
  if( i >= _memory.size() )
  {
    _error_string("getMemory", "", "too small memory");
  }
  return _memory[i];
}

void NonlinearSolver::memory()
{
  _memory.set_size( getNVectors() );
  std::string type = getVisitor()->getVectorType();
  // int level = getVisitor()->getVectorLevel();
  // std::cerr << "IterativeSolverWithVisitor::memory() in " << getClassName() << " visitor = "<<getVisitor()->getClassName() << " " << getNVectors()  << " of type " << type << "\n";
  for(int i = 0; i < _memory.size(); i++)
  {
    std::stringstream ss;
    ss<<getClassName()<<"_"<<_id << "_" <<getVisitor()->getClassName()<<"_memory_"<<i;
    _memory[i] = alat::GhostVector( ss.str(), type);
    // std::cerr << "newVector: " << _memory[i] << "\n";
    getVisitor()->newVector(&_memory[i]);
  }
  // assert(0);
}

/*--------------------------------------------------------------------------*/
bool NonlinearSolver::_checkIteration(perulangan::NewtonOutputData& newtonoutputdata) const
{
  perulanganEnums::newtonstatus& status = newtonoutputdata.newton_status;

  std::stringstream infobuffer;
  infobuffer << " "<< std::setiosflags(std::ios::scientific) << std::setprecision(3) << _newtoninternaldata.dxnorm;
  if(_newtoninternaldata.use_linesearch)
  {
    infobuffer << " " << std::setw(2) << _newtoninternaldata.irelax;
  }
  if(_newtoninternaldata.use_regularization)
  {
    infobuffer << " ["<< std::setiosflags(std::ios::scientific) << std::setprecision(2) << _newtoninternaldata.modelfactor <<" "<<_newtoninternaldata.regularization << "]";
  }
  if(_newtoninternaldata.newmatrix)
  {
    infobuffer << " M";
  }
  else
  {
    infobuffer << " -";
  }
  infobuffer <<  " (" <<_newtoninternaldata.lintol << " " << std::setw(3) << _newtoninternaldata.nliniter << ")\n";
  _info._printbuff = infobuffer.str();
  // std::cerr << "_newtoninternaldata.iteration="<<_newtoninternaldata.iteration<< " _newtoninternaldata.residual="<<_newtoninternaldata.residual<<"\n";
  _info.checkIteration(_newtoninternaldata.status_iteration, _newtoninternaldata.residual, _newtoninternaldata.iteration);
  // std::cerr << "after info.checkIteration() _newtoninternaldata.status_iteration="<<perulanganEnums::iterationStatusToString(_newtoninternaldata.status_iteration)<<"\n";
  if(_newtoninternaldata.residual>std::numeric_limits<double>::max())
  {
    assert(0);
  }
  if(_newtoninternaldata.status_iteration == perulanganEnums::IterationStatusConverged)
  {
    status = perulanganEnums::NewtonStatusConverged;
  }
  else if(_newtoninternaldata.status_iteration == perulanganEnums::IterationStatusDiverged)
  {
    status = perulanganEnums::NewtonStatusDiverged;
  }
  else if(_newtoninternaldata.status_iteration == perulanganEnums::IterationStatusRunning)
  {
    status = perulanganEnums::NewtonStatusRunning;
  }
  else{
    assert(0);
  }
  double rho = _info.getLastReductionRate();
  // std::cerr << "NonlinearSolver::_checkIteration() rho="<<rho<<" _newtoninputdata.rhomatrix="<<_newtoninputdata.rhomatrix<<"\n";
  if(rho < _newtoninputdata.rhomatrix)
  {
    _newtoninternaldata.newmatrix = false;
  }
  else
  {
    _newtoninternaldata.newmatrix = true;
  }

  // if(_newtoninternaldata.iteration)
  // {
  //   double rho = info.getLastReductionRate();
    _newtoninternaldata.lintol = fmin(rho, 1.0)*_newtoninputdata.lineartoleranceincrease;
  // }
  // else
  // {
  //   _lintol = _lineartoleranceincrease;
  // }
  _newtoninternaldata.lintol = fmax( _newtoninternaldata.lintol, 0.99*_info.getMissingRate() );

  // std::cerr << "_newtoninternaldata="<<_newtoninternaldata << "\n";
  // std::cerr << "NonlinearSolver::_checkIteration() status="<<perulanganEnums::newtonStatusToString(status)<<"\n";

  newtonoutputdata.niter = _newtoninternaldata.iteration;

  if( status == perulanganEnums::NewtonStatusConverged or status == perulanganEnums::NewtonStatusDiverged or (_newtoninternaldata.use_linesearch and status == perulanganEnums::NewtonStatusMaxLineSearchAttained) )
  {
    return true;
  }
  else
  {
    if(_newtoninternaldata.newmatrix)
    {
     newtonoutputdata.nredo_matrix++;
    }
   return false;
  }
}

/*--------------------------------------------------------------------------*/
void NonlinearSolver::_lineSearchMonotonicity(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostVector& u, const alat::GhostVector& du, alat::GhostVector& r, const alat::GhostVector& f) const
{
  _newtoninternaldata.irelax = 0;
  _newtoninternaldata.omega = 1.0;
  double omegaold;
  while(_newtoninternaldata.residual > _newtoninputdata.monotonyfactor*_newtoninternaldata.residual_old)
  {
    if(_newtoninternaldata.irelax == _newtoninputdata.maxnlinesearch)
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusMaxLineSearchAttained;
      return;
      // _error_string("_lineSearchMonotonicity", "maxiter linesearch attained", _maxnlinesearch);
    }
    _newtoninternaldata.irelax++;
    if(_newtoninputdata.printlinesearch)
    {
      std::cerr << "irelax="<<_newtoninternaldata.irelax << " res="<<_newtoninternaldata.residual << "resold="<<_newtoninternaldata.residual_old  << "_monotonyfactor*resold="<<_newtoninputdata.monotonyfactor*_newtoninternaldata.residual_old <<"\n";
    }
    omegaold = _newtoninternaldata.omega;
    _newtoninternaldata.omega = omegaold*_newtoninputdata.omegalinesearch;
    // std::cerr << "omegaold="<<omegaold << " omega="<<omega<<"\n";
    getVisitor()->vectorAdd(u, _newtoninternaldata.omega-omegaold, du);
    getVisitor()->residual(newtonoutputdata.residual_status, r, u, f);
    _newtoninternaldata.residual =  sqrt(getVisitor()->vectorDot(r, r));
  }
}

/*--------------------------------------------------------------------------*/
void NonlinearSolver::_lineSearchArmijo(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostVector& u, const alat::GhostVector& du, alat::GhostVector& r, const alat::GhostVector& f, double kappa) const
{
  //
  // on veut  \phi(x+\omega*dx) \le \phi(x) + c1 \omega* dx*\phi'(x)dx
  // avec \phi(x) = \fra12\|b-f(x)\|^2, donc  \phi'(x) = -r*f'(x)dx
  // cela donne
  // res^2 \le resold^2 - c1*\omega\kappa \qquad  \kappa=r*f'(x)dx
  //
  _newtoninternaldata.irelax = 0;
  _newtoninternaldata.omega = 1.0;
  double omegaold;
  bool notok = ( _newtoninternaldata.residual*_newtoninternaldata.residual > _newtoninternaldata.residual_old*_newtoninternaldata.residual_old - 0.0001*_newtoninternaldata.omega*kappa );
  while(notok)
  {
    // std::cerr << "_lineSearchArmijo() res="<<res <<" resold="<<resold<<" kappa="<<kappa<< " omega="<<omega<<"\n";
    if(_newtoninternaldata.irelax == _newtoninputdata.maxnlinesearch)
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusMaxLineSearchAttained;
      return;
      // _error_string("_lineSearchArmijo", "maxiter linesearch attained", _maxnlinesearch);
    }
    _newtoninternaldata.irelax++;
    omegaold = _newtoninternaldata.omega;
    _newtoninternaldata.omega = omegaold*_newtoninputdata.omegalinesearch;
    getVisitor()->vectorAdd(u, _newtoninternaldata.omega-omegaold, du);
    getVisitor()->residual(newtonoutputdata.residual_status, r, u, f);
    _newtoninternaldata.residual =  sqrt(getVisitor()->vectorDot(r, r));
    notok = ( _newtoninternaldata.residual*_newtoninternaldata.residual > _newtoninternaldata.residual_old*_newtoninternaldata.residual_old - 0.0001*_newtoninternaldata.omega*kappa );
  }
}

/*--------------------------------------------------------------------------*/
void NonlinearSolver::_lineSearchWolfe(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostVector& u, const alat::GhostVector& du, alat::GhostVector& r, const alat::GhostVector& f, double kappa) const
{
  // de plus (par rapport à Armijo) on demande
  //    \phi'(x+omega*dx)dx \ge c2   \phi'(x)dx
  //    \|\phi'(x+omega*dx)dx\| \le c2   \|\phi'(x)dx\|
  _newtoninternaldata.irelax = 0;
  _newtoninternaldata.omega = 1.0;
  double omegaold;
  alat::GhostVector& h = getMemory(2);
  perulanganEnums::residualstatus status_linearization = perulanganEnums::ResidualStatusOk;

  getVisitor()->vectorZero(h);
  getVisitor()->computeLinearization(newtonoutputdata.residual_status, h, u, du);
  if(newtonoutputdata.residual_status != perulanganEnums::ResidualStatusOk)
  {
    _error_string( "_lineSearchWolfe", "linearization not ok", perulanganEnums::residualStatusToString(status_linearization) );
  }
  double kappanew = getVisitor()->vectorDot(r, h);
  // std::cerr << "_lineSearchArmijo() res="<<res <<" resold="<<resold<<" kappa="<<kappa<<" kappanew="<<kappanew<< " omega="<<omega<<"\n";

  // bool notok = ( res*res > resold*resold - 0.0001*omega*kappa ) or ( kappanew < 0.1*kappa );
  bool notok = ( _newtoninternaldata.residual*_newtoninternaldata.residual > _newtoninternaldata.residual_old*_newtoninternaldata.residual_old - 0.0001*_newtoninternaldata.omega*kappa )or ( fabs(kappanew) < 0.1*fabs(kappa) );
  while(notok)
  {
    // std::cerr << "_lineSearchArmijo() res="<<res <<" resold="<<resold<<" kappa="<<kappa<<" kappanew="<<kappanew<< " omega="<<omega<<"\n";
    if(_newtoninternaldata.irelax == _newtoninputdata.maxnlinesearch)
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusMaxLineSearchAttained;
      return;
      // _error_string("_lineSearchWolfe", "maxiter linesearch attained", _maxnlinesearch);
    }
    _newtoninternaldata.irelax++;
    omegaold = _newtoninternaldata.omega;
    _newtoninternaldata.omega = omegaold*_newtoninputdata.omegalinesearch;
    getVisitor()->vectorAdd(u, _newtoninternaldata.omega-omegaold, du);
    getVisitor()->residual(newtonoutputdata.residual_status, r, u, f);
    _newtoninternaldata.residual =  sqrt(getVisitor()->vectorDot(r, r));
    getVisitor()->vectorZero(h);
    getVisitor()->computeLinearization(newtonoutputdata.residual_status, h, u, du);
    if(newtonoutputdata.residual_status != perulanganEnums::ResidualStatusOk)
    {
      newtonoutputdata.newton_status = perulanganEnums::NewtonStatusLinearNotOk;
      assert(0);
      return;
      // _error_string( "_lineSearchWolfe", "linearization not ok", perulanganEnums::residualStatusToString(status_linearization) );
    }
    kappanew = getVisitor()->vectorDot(r, h);
    notok = ( _newtoninternaldata.residual*_newtoninternaldata.residual > _newtoninternaldata.residual_old*_newtoninternaldata.residual_old - 0.0001*_newtoninternaldata.omega*kappa )or ( fabs(kappanew) < 0.1*fabs(kappa) );
  }
}

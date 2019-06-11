#include  "Solvers/exactsolutions.hpp"
#include  "application.hpp"
#include  "model.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
RDcosh::RDcosh(double diff, double alpha) : solvers::FunctionInterface()
{
  _b = sqrt(alpha/diff);
  _c = -1.0/cosh(_b);
}
std::string RDcosh::getClassName() const {return "RDcosh";}
void RDcosh::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = 1.0+_c*cosh(_b*x);
}
void RDcosh::x(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _b*_c*sinh(_b*x);
}
void RDcosh::xx(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _b*_b*_c*cosh(_b*x);
}
void RDcosh::y (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void RDcosh::z (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void RDcosh::t (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void RDcosh::yy(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void RDcosh::zz(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}

/*--------------------------------------------------------------------------*/
CDExpLayer::CDExpLayer(double diff) : solvers::FunctionInterface()
{
  _diff = diff;
}
std::string CDExpLayer::getClassName() const {return "CDExpLayer";}
void CDExpLayer::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = (1.0+y)*( 1.0-exp((x-1)/_diff) )/( 1.0-exp(-2.0/_diff) );
  u[0] = (1.0-exp((1.0-x)/_diff))/(1.0-exp(2.0/_diff));
}
void CDExpLayer::x(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = -(1.0+y)*exp((x-1)/_diff)/( 1.0-exp(-2.0/_diff) )/_diff;
  u[0] = exp((1.0-x)/_diff)/(1.0-exp(2.0/_diff))/_diff;
}
void CDExpLayer::y(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = ( 1.0-exp((x-1)/_diff) )/( 1.0-exp(-2.0/_diff) );
  u[0] = ( exp(-2.0/_diff)-exp((-1.0-x)/_diff) )/( exp(-2.0/_diff)-1.0 );
}
void CDExpLayer::xx(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void CDExpLayer::z (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void CDExpLayer::t (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void CDExpLayer::yy(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void CDExpLayer::zz(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}

/*--------------------------------------------------------------------------*/
Application::~Application() {}
Application::Application(std::string exactsolution, bool strongdir): solvers::ApplicationInterface(), _exactsolution(exactsolution), _localmodel(NULL), _strongdir(strongdir){}
Application::Application( const Application& application): solvers::ApplicationInterface(application)
{
  assert(0);
}
Application& Application::operator=( const Application& application)
{
  assert(0);
  solvers::ApplicationInterface::operator=(application);
  return *this;
}
std::string Application::getClassName() const
{
  return "Application";
}
/*--------------------------------------------------------------------------*/
void Application::initApplication(const mesh::MeshUnitInterface* mesh, const alat::StringVector& varnames, const alat::StringVector& varnamesdata, const alat::armaivec& ncomps, const alat::armaivec& ncompsdata)
{
  _localmodel = dynamic_cast<const Model*>(_model);
  assert(_localmodel);
  ApplicationInterface::initApplication(mesh, varnames, varnamesdata, ncomps, ncompsdata);
}
/*--------------------------------------------------------------------------*/
std::shared_ptr<solvers::FunctionInterface> Application::newExactSolution(std::string varname) const
{
  if(_exactsolution=="none") return solvers::ApplicationInterface::newExactSolution(varname);
  if(_exactsolution=="constant")
  {
    double a=7.0;
    return std::unique_ptr<solvers::FunctionInterface>(new solvers::ConstantSolution(a));
  }
  else if(_exactsolution=="linear")
  {
    double a=7.0, b = 1.2, c = 1.7, d = 1.9;
    if(_mesh->getDimension()<3) {d=0.0;}
    if(_mesh->getDimension()<2) {c=0.0;}
    return std::unique_ptr<solvers::FunctionInterface>(new solvers::LinearSolution(a, b, c, d));
  }
  else if(_exactsolution=="quadratic")
  {
    // double a=2.0, b = 1.2, c = 1.7, d = 1.9, e = 0.1, f = 0.2, g=0.3, h=0.4, i=0.5;
    double a=_mesh->getDimension(), b = 0.0, c = 0.0, d = -1.0, e = -1.0, f = 0.0, g=0.0, h=0.0, i=-1.0;
    if(_mesh->getDimension()<3) {g=h=i=0.0;}
    if(_mesh->getDimension()<2) {c=e=f=0.0;}
    return std::unique_ptr<solvers::FunctionInterface>(new solvers::QuadraticSolution(a, b, c, d, e, f, g, h, i));
  }
  else if(_exactsolution=="cosinus")
  {
    double a=1.0, b = 1.2, c = 1.7, d = 1.9;
    if(_mesh->getDimension()<3) {d=0.0;}
    if(_mesh->getDimension()<2) {c=0.0;}
    return std::unique_ptr<solvers::FunctionInterface>(new solvers::CosinusSolution(a, b, c, d));
  }
  else if(_exactsolution=="rdcosh")
  {
    double diff=_localmodel->_diff;
    double alpha=_localmodel->_alpha;
    return std::unique_ptr<solvers::FunctionInterface>(new RDcosh(diff, alpha));
  }
  else if(_exactsolution=="cdexplayer")
  {
    double diff=_localmodel->_diff;
    assert(_localmodel->_alpha==0.0);
    return std::unique_ptr<solvers::FunctionInterface>(new CDExpLayer(diff));
  }
  else
  {
    _error_string("Application","unknwon exactsolution",_exactsolution);
  }
}
/*--------------------------------------------------------------------------*/
bool Application::isStrongDirichlet(int color)const
{
  if(not _strongdir) {return false;}
  if(_exactsolution=="cdexplayer")
  {
    if( color==111 or color==333)
    {
      return false;
    }
  }
  return true;
}
class RightHandSideExactSolution : public solvers::RightHandSideInterface{
private:
  const Model* localmodel;
  const solvers::FunctionInterface& solution;
public:
  RightHandSideExactSolution(const Model* localmode_in, const solvers::FunctionInterface& solution_in) : solvers::RightHandSideInterface(), localmodel(localmode_in), solution(solution_in){}
  std::string getClassName()const {return "RightHandSideExactSolution";}
  void operator()(alat::armavec& f, double x, double y, double z, double t)const
  {
    alat::armavec u(1), ux(1), uy(1), uz(1), uxx(1), uyy(1), uzz(1);
    arma::mat r(1,1), um(1,1);
    solution(u, x, y, z);
    solution.x(ux, x, y, z);
    solution.y(uy, x, y, z);
    solution.z(uz, x, y, z);
    solution.xx(uxx, x, y, z);
    solution.yy(uyy, x, y, z);
    solution.zz(uzz, x, y, z);
    um(0,0) = u[0];
    // u[0] = solution[0]->operator()(x,y,z);
    localmodel->reaction(r.col(0), um.col(0));
    double d = localmodel->diffusion(x,y,z);
    alat::armavec beta(3);
    localmodel->beta(beta, x, y, z);
    f[0] = r(0,0);
    // f[0] += beta[0]*solution[0]->x(x,y,z)+beta[1]*solution[0]->y(x,y,z)+beta[1]*solution[0]->z(x,y,z);
    // f[0] -= d*(solution[0]->xx(x,y,z)+solution[0]->yy(x,y,z)+solution[0]->zz(x,y,z));
    f[0] += beta[0]*ux[0]+beta[1]*uy[0]+beta[1]*uz[0];
    f[0] -= d*(uxx[0]+uyy[0]+uzz[0]);
  }
};
class NeumannExactSolution : public solvers::NeumannInterface{
private:
  const Model* localmodel;
  const solvers::FunctionInterface& solution;
public:
  NeumannExactSolution(const Model* localmode_in, const solvers::FunctionInterface& solution_in) : solvers::NeumannInterface(), localmodel(localmode_in), solution(solution_in){}
  std::string getClassName()const {return "NeumannExactSolution";}
  void operator()(alat::armavec& f, double nx, double ny, double nz, double x, double y, double z, double t)const
  {
    alat::armavec ux(1), uy(1), uz(1);
    arma::mat r(1,1), um(1,1);
    solution.x(ux, x, y, z);
    solution.y(uy, x, y, z);
    solution.z(uz, x, y, z);
    double d = localmodel->diffusion(x,y,z);
    f[0] = d*(nx*ux[0]+ny*uy[0]+nz*uz[0]);
  }
};

/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::RightHandSideInterface> Application::newRightHandSide(std::string varname) const
{
  assert(_exactsolutions.size()==1);
  assert(_localmodel);
  if( _exactsolution=="cdexplayer")
  {
    return std::unique_ptr<solvers::RightHandSideInterface>(new solvers::RightHandSideConstant(0.0));
  }
  else if( _exactsolution=="rdcosh")
  {
    return std::unique_ptr<solvers::RightHandSideInterface>(new solvers::RightHandSideConstant(1.0));
  }
  return std::unique_ptr<solvers::RightHandSideInterface>(new RightHandSideExactSolution(_localmodel, getExactSolution(0)));
}
/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::DirichletInterface> Application::newDirichlet(std::string varname) const
{
  assert(_exactsolutions.size()==1);
  if( _exactsolution=="rdcosh")
  {
    return std::unique_ptr<solvers::DirichletInterface>(new solvers::DirichletConstant(0.0));
  }
  return std::unique_ptr<solvers::DirichletInterface>(new solvers::DirichletExactSolution(getExactSolution(0)));
}
/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::NeumannInterface> Application::newNeumann(std::string varname) const
{
  assert(_exactsolutions.size()==1);
  assert(_localmodel);
  return std::unique_ptr<solvers::NeumannInterface>(new NeumannExactSolution(_localmodel, getExactSolution(0)));
}
/*--------------------------------------------------------------------------*/
std::shared_ptr<solvers::FunctionInterface> Application::newDataFunction(std::string varname) const
{
  assert(_datafunctions.size()==1);
  return _localmodel->getBeta();
}

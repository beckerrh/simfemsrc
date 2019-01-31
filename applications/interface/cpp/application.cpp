#include  "Solvers/exactsolutions.hpp"
#include  "application.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
class InterfacePositionStraight : public solvers::FunctionInterface
{
protected:
  double _xgamma;
public:
  InterfacePositionStraight(double xgamma) : solvers::FunctionInterface(), _xgamma(xgamma) {}
  void operator()(alat::armavec& u, double x, double y, double z, double t)const
  {
    u[0] = x- _xgamma;
  }
};
/*--------------------------------------------------------------------------*/
class InterfacePositionCircle : public solvers::FunctionInterface
{
public:
  void operator()(alat::armavec& u, double x, double y, double z, double t)const
  {
    u[0] = x*x + y*y + z*z - 0.5;
  }
};

/*--------------------------------------------------------------------------*/
CutExactSolutionQuadraticCircle::CutExactSolutionQuadraticCircle(double k1, double k2) : solvers::FunctionInterface(), _k1(k1), _k2(k2), _phi(NULL){}
void CutExactSolutionQuadraticCircle::setPhi(const solvers::FunctionInterface* phi)
{
  _phi = phi;
  _p.set_size(1);
  (*phi)(_p, 0.0, 0.0, 0.0);
  _r02 = -_p[0];
}

std::string CutExactSolutionQuadraticCircle::getName() const {return "CutExactSolutionQuadraticCircle";}
void CutExactSolutionQuadraticCircle::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  (*_phi)(_p, x, y, z);
  if(_p[0]<0.0)
  {
    u[0] = (x*x+y*y)/_k1;
  }
  else
  {
    u[0] = (x*x+y*y)/_k2 - _r02/_k2 + _r02/_k1;
  }
}
void CutExactSolutionQuadraticCircle::x (alat::armavec& u, double x, double y, double z, double t) const
{
  (*_phi)(_p, x, y, z);
  if(_p[0]<0.0)
  {
    u[0] = 2.0*x/_k1;
  }
  else
  {
    u[0] = 2.0*x/_k2;
  }
}
void CutExactSolutionQuadraticCircle::y (alat::armavec& u, double x, double y, double z, double t) const
{
  (*_phi)(_p, x, y, z);
  if(_p[0]<0.0)
  {
    u[0] = 2.0*y/_k1;
  }
  else
  {
    u[0] = 2.0*y/_k2;
  }
}
void CutExactSolutionQuadraticCircle::z (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void CutExactSolutionQuadraticCircle::t (alat::armavec& u, double x, double y, double z, double t) const{assert(0);}
void CutExactSolutionQuadraticCircle::xx(alat::armavec& u, double x, double y, double z, double t) const{assert(0);}
void CutExactSolutionQuadraticCircle::yy(alat::armavec& u, double x, double y, double z, double t) const{assert(0);}
void CutExactSolutionQuadraticCircle::zz(alat::armavec& u, double x, double y, double z, double t) const{assert(0);}

/*--------------------------------------------------------------------------*/
CutExactSolutionLinearStraight::CutExactSolutionLinearStraight(double xgamma, double k1, double k2) : solvers::FunctionInterface(), _xgamma(xgamma), _k1(k1), _k2(k2)
{
  _p1 = 2.0*k2/( (k2-k1)*xgamma + k1+k2);
  _p2 = 2.0*k1/( (k2-k1)*xgamma + k1+k2);
}

std::string CutExactSolutionLinearStraight::getName() const {return "CutExactSolutionLinearStraight";}
void CutExactSolutionLinearStraight::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  if(x<=_xgamma)
  {
    u[0] = -1.0 + _p1*(x+1.0) + 0.0*y;
  }
  else
  {
    u[0] = 1.0 + _p2*(x-1.0) + 0.0*y;
  }
}
void CutExactSolutionLinearStraight::x (alat::armavec& u, double x, double y, double z, double t) const
{
  if(x<=_xgamma)
  {
    u[0] = _p1;
  }
  else
  {
    u[0] = _p2;
  }
}
void CutExactSolutionLinearStraight::y (alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = 0.0;
}
void CutExactSolutionLinearStraight::z (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void CutExactSolutionLinearStraight::t (alat::armavec& u, double x, double y, double z, double t) const{assert(0);}
void CutExactSolutionLinearStraight::xx(alat::armavec& u, double x, double y, double z, double t) const{assert(0);}
void CutExactSolutionLinearStraight::yy(alat::armavec& u, double x, double y, double z, double t) const{assert(0);}
void CutExactSolutionLinearStraight::zz(alat::armavec& u, double x, double y, double z, double t) const{assert(0);}


/*--------------------------------------------------------------------------*/
CutExactSolutionQuadraticStraight::CutExactSolutionQuadraticStraight(double xgamma, double k1, double k2) : solvers::FunctionInterface(), _xgamma(xgamma), _k1(k1), _k2(k2)
{
}

std::string CutExactSolutionQuadraticStraight::getName() const {return "CutExactSolutionQuadraticStraight";}
void CutExactSolutionQuadraticStraight::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  if(x<=_xgamma)
  {
    u[0] = x*x/_k1;
  }
  else
  {
    u[0] = (x*x-_xgamma*_xgamma)/_k2 + _xgamma*_xgamma/_k1;
  }
}
void CutExactSolutionQuadraticStraight::x (alat::armavec& u, double x, double y, double z, double t) const
{
  if(x<=_xgamma)
  {
    u[0] = 2.0*x/_k1;
  }
  else
  {
    u[0] = 2.0*x/_k2;
  }
}
void CutExactSolutionQuadraticStraight::y (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void CutExactSolutionQuadraticStraight::z (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void CutExactSolutionQuadraticStraight::t (alat::armavec& u, double x, double y, double z, double t) const{assert(0);}
void CutExactSolutionQuadraticStraight::xx(alat::armavec& u, double x, double y, double z, double t) const{assert(0);}
void CutExactSolutionQuadraticStraight::yy(alat::armavec& u, double x, double y, double z, double t) const{assert(0);}
void CutExactSolutionQuadraticStraight::zz(alat::armavec& u, double x, double y, double z, double t) const{assert(0);}


/*--------------------------------------------------------------------------*/
Application::~Application() {}
Application::Application(std::string applicationname): solvers::ApplicationInterface(), _localmodel(NULL), _applicationname(applicationname){}
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

  if(hasExactSolution("U"))
  {
    solvers::FunctionInterface& uex = getExactSolution("U");
    CutExactSolutionQuadraticCircle* cutexactsolutionquadraticcircle = dynamic_cast<CutExactSolutionQuadraticCircle*>(&uex);
    if(cutexactsolutionquadraticcircle)
    {
      const solvers::FunctionInterface& phi = getDataFunction("Phi");
      cutexactsolutionquadraticcircle->setPhi(&phi);
    }    
  }
}

/*--------------------------------------------------------------------------*/
bool Application::isStrongDirichlet(int color)const
{
  return true;
}
/*--------------------------------------------------------------------------*/
std::shared_ptr<solvers::FunctionInterface> Application::newDataFunction(std::string varname) const
{
  // std::cerr << "_applicationname="<<_localmodel->_xgamma<<"\n";
  if(_applicationname=="linear_straight" or _applicationname=="quadratic_straight")
  {
    return std::unique_ptr<solvers::FunctionInterface>(new InterfacePositionStraight(_localmodel->_xgamma));
  }
  else if(_applicationname=="quadratic_straight")
  {
    return std::unique_ptr<solvers::FunctionInterface>(new InterfacePositionCircle);
  }
  return std::unique_ptr<solvers::FunctionInterface>(nullptr);
}
/*--------------------------------------------------------------------------*/
std::shared_ptr<solvers::FunctionInterface> Application::newExactSolution(std::string varname) const
{
  if(_applicationname=="constant")
  {
    double a=7.0;
    return std::unique_ptr<solvers::FunctionInterface>(new solvers::ConstantSolution(a));
  }
  else if(_applicationname=="linear")
  {
    double a=7.0, b = 1.0, c = 2.0, d = 1.0;
    if(_mesh->getDimension()<3) {d=0.0;}
    if(_mesh->getDimension()<2) {c=0.0;}
    return std::unique_ptr<solvers::FunctionInterface>(new solvers::LinearSolution(a, b, c, d));
  }
  else if(_applicationname=="linear_straight")
  {
    assert(_localmodel);
    return std::unique_ptr<solvers::FunctionInterface>(new CutExactSolutionLinearStraight(_localmodel->_xgamma, _localmodel->_kin, _localmodel->_kex));
  }
  else if(_applicationname=="quadratic_straight")
  {
    assert(_localmodel);
    return std::unique_ptr<solvers::FunctionInterface>(new CutExactSolutionQuadraticStraight(_localmodel->_xgamma, _localmodel->_kin, _localmodel->_kex));
  }
  else if(_applicationname=="quadratic_circle")
  {
    assert(_localmodel);
    return std::unique_ptr<solvers::FunctionInterface>(new CutExactSolutionQuadraticCircle(_localmodel->_kin, _localmodel->_kex));
  }
  else
  {
    return std::unique_ptr<solvers::FunctionInterface>(nullptr);
  }
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
    assert(0);
    double d = localmodel->_kin;
    alat::armavec uxx(1), uyy(1), uzz(1);
    solution.xx(uxx, x, y, z);
    solution.yy(uyy, x, y, z);
    solution.zz(uzz, x, y, z);
    f[0] -= d*(uxx[0]+uyy[0]+uzz[0]);
  }
};

/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::RightHandSideInterface> Application::newRightHandSide(std::string varname) const
{
  assert(_localmodel);
  if(_applicationname=="cardioide")
  {
    return std::unique_ptr<solvers::RightHandSideInterface>(new solvers::RightHandSideConstant(1));
  }

  if( _applicationname=="constant" or  _applicationname=="linear" or  _applicationname=="linear_straight")
  {
    return std::unique_ptr<solvers::RightHandSideInterface>(new solvers::RightHandSideConstant());
  }
  else if(_applicationname=="quadratic_straight")
  {
    return std::unique_ptr<solvers::RightHandSideInterface>(new solvers::RightHandSideConstant(-2));
  }
  else if( _applicationname=="quadratic_circle")
  {
    return std::unique_ptr<solvers::RightHandSideInterface>(new solvers::RightHandSideConstant(-4));
  }
  return std::unique_ptr<solvers::RightHandSideInterface>(new RightHandSideExactSolution(_localmodel,getExactSolution(0)));
}
std::unique_ptr<solvers::DirichletInterface> Application::newDirichlet(std::string varname) const
{
  if(_applicationname=="cardioide")
  {
    return std::unique_ptr<solvers::DirichletInterface>(new solvers::DirichletConstant());
  }
  return std::unique_ptr<solvers::DirichletInterface>(new solvers::DirichletExactSolution(getExactSolution(0)));
}

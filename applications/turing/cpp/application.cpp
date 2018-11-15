#include  "application.hpp"
#include  "model.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
class InitialCondition : public solvers::InitialConditionInterface
{
private:
  double _a, _b;

public:
  InitialCondition(double a, double b) : solvers::InitialConditionInterface(){
    _a = a; _b = b;
  }
  void operator()(arma::vec& u, double x, double y, double z, double t)const{
    u.fill(arma::fill::zeros);
    u[0] = _a;
    u[1] = _b/_a;
    double d = 1.1;
    if(_dim==1)
    {
      if( (x>=-0.4) and (x<=0.0))
      {
        u[0] *= d;
      }
      if( (x>=-0.2) and (x<=0.2))
      {
        u[1] *= d;
      }
      return;
    }
    if(_dim==2)
    {
      if( (x>=-0.4) and (x<=0.0) and (y>=-0.4) and (y<=0.0) )
      {
        u[0] *= d;
      }
      if( (x>=-0.2) and (x<=0.2) and (y>=-0.2) and (y<=0.2) )
      {
        u[1] *= d;
      }
      return;
    }
    if( (x>=-0.4) and (x<=0.0) and (y>=-0.4) and (y<=0.0) and (z>=-0.4) and (z<=0.0) )
    {
      u[0] *= d;
    }
    if( (x>=-0.2) and (x<=0.2) and (y>=-0.2) and (y<=0.2) and (z>=-0.2) and (z<=0.2) )
    {
      u[1] *= d;
    }
  }
};

/*--------------------------------------------------------------------------*/
Application::~Application() {}
Application::Application(): solvers::ApplicationInterface(){}
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
std::unique_ptr<solvers::InitialConditionInterface> Application::newInitialCondition(std::string varname) const
{
  const Model* model = dynamic_cast<const Model*>(_model);
  assert(model);
  return std::unique_ptr<solvers::InitialConditionInterface>(new InitialCondition(model->a, model->b));
}

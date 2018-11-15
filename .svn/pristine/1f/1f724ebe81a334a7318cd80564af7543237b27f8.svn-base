#include  "Solvers/application.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
Application::~Application() {}
Application::Application(): ApplicationInterface(){}
Application::Application( const Application& application): ApplicationInterface(application)
{
assert(0);
}
Application& Application::operator=( const Application& application)
{
  assert(0);
  ApplicationInterface::operator=(application);
  return *this;
}
std::string Application::getClassName() const
{
  return "Application";
}

/*--------------------------------------------------------------------------*/
bool Application::isStrongDirichlet(int color)const {return true;}

class Dirichlet : public DirichletInterface
{
  void operator()(arma::vec& u, double x, double y, double z, double t)const{
    u[0] = 3.0*x + 2.0*y + z;
  }
};
std::unique_ptr<solvers::DirichletInterface> Application::newDirichlet(std::string varname) const
{
  return std::unique_ptr<solvers::DirichletInterface>(new Dirichlet);
}

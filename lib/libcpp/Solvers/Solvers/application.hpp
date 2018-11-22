#ifndef __Solvers_Application_hpp
#define __Solvers_Application_hpp

#include  "applicationinterface.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
namespace solvers
{
  class Application : public ApplicationInterface
  {
  protected:

  public:
    ~Application();
    Application();
    Application( const Application& application);
    Application& operator=( const Application& application);
    std::string getClassName() const;

    bool isStrongDirichlet(int color)const;

    std::unique_ptr<solvers::DirichletInterface> newDirichlet(std::string varname) const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

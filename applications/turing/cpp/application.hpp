#ifndef ___Application_hpp
#define ___Application_hpp

#include  "Solvers/applicationinterface.hpp"

/*--------------------------------------------------------------------------*/
  class Application : public solvers::ApplicationInterface
  {
  public:
    ~Application();
    Application();
    Application( const Application& application);
    Application& operator=( const Application& application);
    std::string getClassName() const;

    std::unique_ptr<solvers::InitialConditionInterface> newInitialCondition(std::string varname) const;
    bool isStrongDirichlet(int color)const {return false;}
  };

/*--------------------------------------------------------------------------*/
#endif

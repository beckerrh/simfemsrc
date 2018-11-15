#ifndef ___Solver_hpp
#define ___Solver_hpp

#include  "Solvers/solver.hpp"

/*--------------------------------------------------------------------------*/
  class Solver : public solvers::Solver
  {
  protected:
  public:
    ~Solver();
    Solver();
    Solver( const Solver& solver);
    Solver& operator=( const Solver& solver);
    std::string getClassName() const;

    std::unique_ptr<solvers::FemInterface> newFem(const std::string& varname, solverEnums::fem::femtype fem) const;
    std::unique_ptr<solvers::PdePartInterface> newPdePart(std::string pdepartname) const;
    std::unique_ptr<solvers::ApplicationInterface> newApplication() const;
    std::unique_ptr<solvers::ModelInterface> newModel() const;
    void defineVariablesAndPdeParts();
    solvers::ErrorsMap run();
  };

/*--------------------------------------------------------------------------*/
#endif

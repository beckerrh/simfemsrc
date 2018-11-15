#ifndef ___Solver_hpp
#define ___Solver_hpp

#include  "Solvers/solver.hpp"

/*--------------------------------------------------------------------------*/
  class Solver : public solvers::Solver
  {
  public:
    ~Solver();
    Solver();
    Solver( const Solver& solver);
    Solver(std::shared_ptr<mesh::MeshInterface> mesh=nullptr, solver_options::opts opts=solver_options::none);
    Solver& operator=( const Solver& solver);
    std::string getClassName() const;

    // std::unique_ptr<solvers::PdePartInterface> newPdePart(const solvers::Variable& var) const;
    // void init();
    std::unique_ptr<solvers::ApplicationInterface> newApplication() const;
    std::unique_ptr<solvers::ModelInterface> newModel() const;
    void defineVariablesAndPdeParts();
  };

/*--------------------------------------------------------------------------*/
#endif

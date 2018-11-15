#ifndef ___Solver_hpp
#define ___Solver_hpp

#include  "Solvers/solver.hpp"

/*--------------------------------------------------------------------------*/
  class Solver : public solvers::Solver
  {
  protected:
    std::unique_ptr<solvers::MeshUnitWithDataInterface>newMeshUnitWithDataPlain(const mesh::MeshUnitInterface* mesh) const;
    void initParameters();

  public:
    ~Solver();
    Solver( const Solver& solver);
    Solver(std::shared_ptr<mesh::MeshInterface> mesh=nullptr, solver_options::opts opts=solver_options::none);
    Solver& operator=( const Solver& solver);
    std::string getClassName() const;

    std::unique_ptr<solvers::FemInterface> newFem(const std::string& varname, solverEnums::fem::femtype fem) const;
    std::unique_ptr<solvers::ApplicationInterface> newApplication() const;
    std::unique_ptr<solvers::ModelInterface> newModel() const;
    void defineVariablesAndPdeParts();
    solvers::ErrorsMap run();
    void init();
  };

/*--------------------------------------------------------------------------*/
#endif

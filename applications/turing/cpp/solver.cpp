#include  "solver.hpp"
#include  "application.hpp"
#include  "pdepart.hpp"
#include  "model.hpp"
#include  "Solvers/p1.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
Solver::~Solver() {}
Solver::Solver(): solvers::Solver()
{
  // parameters.doubles["b"] = 2.1;
  _parameters.doubles["a"] = 1.0;
  _parameters.doubles["b"] = 1.9;
  _parameters.doubles["k0"] = 0.002;
  _parameters.doubles["k1"] = 0.2;
}
Solver::Solver( const Solver& solver): solvers::Solver(solver)
{
  assert(0);
}
Solver::Solver(std::shared_ptr<mesh::MeshInterface> mesh, solver_options::opts opts) : solvers::Solver(mesh, opts)
{
  _parameters.doubles["a"] = 1.0;
  _parameters.doubles["b"] = 1.9;
  _parameters.doubles["k0"] = 0.002;
  _parameters.doubles["k1"] = 0.2;
}
Solver& Solver::operator=( const Solver& solver)
{
  assert(0);
  solvers::Solver::operator=(solver);
  return *this;
}
std::string Solver::getClassName() const
{
  return "Solver";
}
/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::ApplicationInterface> Solver::newApplication() const
{
  return std::unique_ptr<solvers::ApplicationInterface>(new Application());

}
std::unique_ptr<solvers::ModelInterface> Solver::newModel() const
{
  return std::unique_ptr<solvers::ModelInterface>(new Model());
}
/*--------------------------------------------------------------------------*/
void Solver::defineVariablesAndPdeParts()
{
  addVariablePlain("U", 2, solverEnums::fem::P1);

  std::unique_ptr<solvers::PdePartInterface> pdepart;
  pdepart = std::unique_ptr<solvers::PdePartInterface>(new PdePart("U"));
  addPdePartPlain("laplace", std::move(pdepart));

  perulangan::NewtonInputData newtoninputdata;
  // newtoninputdata.rhomatrix = 0.1;
  getNonlinearSolver()->setNewtonInputData(newtoninputdata);

  solvers::solverdata solverdata;
  solverdata.timestepratio=0.01;
  // solverdata.rhomatrix=0.1;
  setSolverData(solverdata);
}

#include  "solver.hpp"
#include  "application.hpp"
#include  "model.hpp"
#include  "traditionalintegration.hpp"
#include  "nitscheintegration.hpp"
#include  "newnitscheintegration.hpp"
#include  "traditional.hpp"
#include  "nitsche.hpp"
#include  "newnitsche.hpp"
#include  "Solvers/p1.hpp"
#include  "Solvers/cr1.hpp"
#include  "Solvers/rt0.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
Solver::~Solver() {}
Solver::Solver(): solvers::Solver()
{
  _parameters.strings["method"] = "none";
  _parameters.strings["application"] = "none";
  _parameters.strings["beta"] = "none";
  _parameters.strings["fem"] = "none";
  _parameters.doubles["alpha"] = 0.0;
  _parameters.doubles["diff"] = 1.0;
  _parameters.doubles["deltasupg"] = 0.0;
  _parameters.doubles["gamma"] = 2.0;
  _parameters.bools["lumpedmass"] = true;
  _parameters.bools["symmetric"] = true;
}
Solver::Solver( const Solver& solver): solvers::Solver(solver)
{
assert(0);
}
std::string Solver::getClassName() const
{
  return "Solver";
}
/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::ApplicationInterface> Solver::newApplication() const
{
  std::string method = _parameters.strings["method"];
  bool strongdir=true;
  if(method=="nitsche" or method=="newnitsche" or method=="nitscheintegration" or method=="newnitscheintegration")
  {
    strongdir = false;
  }
  std::string application = _parameters.strings["application"];
  return std::unique_ptr<solvers::ApplicationInterface>(new Application(application, strongdir));
}
std::unique_ptr<solvers::ModelInterface> Solver::newModel() const
{
  return std::unique_ptr<solvers::ModelInterface>(new Model());
}
/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::FemInterface> Solver::newFem(const std::string& varname, solverEnums::fem::femtype fem) const
{
  if(fem == solverEnums::fem::P1)
  {
    return std::unique_ptr<solvers::FemInterface>(new solvers::P1());
  }
  else if(fem == solverEnums::fem::CR1)
  {
    return std::unique_ptr<solvers::FemInterface>(new solvers::CR1());
  }
  return solvers::Solver::newFem(varname, fem);
}

/*--------------------------------------------------------------------------*/
void Solver::defineVariablesAndPdeParts()
{
  std::string method = _parameters.strings["method"];
  std::string fem = _parameters.strings["fem"];

  setOutputOptions(solver_options::solutionfilename) = method;
  if(fem=="P1")
  {
    addVariablePlain("U", 1, solverEnums::fem::P1);
  }
  else if(fem=="CR1")
  {
    addVariablePlain("U", 1, solverEnums::fem::CR1);
  }
  else if(fem=="P2")
  {
    addVariablePlain("U", 1, solverEnums::fem::P2);
  }
  else
  {
    assert(0);
  }
  addDataVariablePlain("beta", 1, solverEnums::fem::RT0);

  addPdePartPlain("laplace", newPdePart(method));

  perulangan::NewtonInputData newtoninputdata;
  // newtoninputdata.rhomatrix = 0.1;
  getNonlinearSolver()->setNewtonInputData(newtoninputdata);

  solvers::solverdata solverdata;
  solverdata.timestepratio=0.01;
  // solverdata.rhomatrix=0.1;
  setSolverData(solverdata);
}

/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::PdePartInterface> Solver::newPdePart(std::string method) const
{
  alat::StringList vars("U");
  std::unique_ptr<solvers::PdePartInterface> pdepart;
  if(method=="traditionalintegration")
  {
    pdepart = std::unique_ptr<solvers::PdePartInterface>(new TraditionalIntegration(vars));
  }
  else if(method=="nitscheintegration")
  {
    pdepart = std::unique_ptr<solvers::PdePartInterface>(new NitscheIntegration(vars));
  }
  else if(method=="newnitscheintegration")
  {
    pdepart = std::unique_ptr<solvers::PdePartInterface>(new NewNitscheIntegration(vars));
  }
  else if(method=="traditional")
  {
    pdepart = std::unique_ptr<solvers::PdePartInterface>(new Traditional(vars));
  }
  else if(method=="nitsche")
  {
    pdepart = std::unique_ptr<solvers::PdePartInterface>(new Nitsche(vars));
  }
  else if(method=="newnitsche")
  {
    pdepart = std::unique_ptr<solvers::PdePartInterface>(new NewNitsche(vars));
  }
  else
  {
    _error_string("newPdePart","unknown bdry type", method);
  }
  return std::move(pdepart);
}
/*--------------------------------------------------------------------------*/
solvers::ErrorsMap Solver::run()
{
  solvers::Solver::run();
  solver_options::errors::opts opts = solver_options::errors::H1;
  opts += solver_options::errors::L2;
  opts += solver_options::errors::L1;
  opts += solver_options::errors::Linf;
  // std::cerr << "ops: " << opts << "\n";
  solvers::ErrorsMap errors = computeErrors("U", opts);
  // std::cerr << "errors: " << errors << "\n";
  return errors;
}

#include  "application.hpp"
#include  "cr1iim.hpp"
#include  "p1cut.hpp"
#include  "p1interface.hpp"
#include  "p1hansbo.hpp"
#include  "p1iim.hpp"
#include  "model.hpp"
#include  "pdepartcutinterface.hpp"
#include  "pdepart.hpp"
#include  "pdepartiim.hpp"
#include  "pdepartmultiplier.hpp"
#include  "solver.hpp"
#include  "meshunitwithdata.hpp"
#include  "Solvers/p1.hpp"
#include  "Mesh/cutinterface.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
Solver::~Solver() {}
// Solver::Solver(): solvers::Solver(){}
Solver::Solver( const Solver& solver): solvers::Solver(solver)
{
  assert(0);
}
Solver::Solver(std::shared_ptr<mesh::MeshInterface> mesh, solver_options::opts opts) : solvers::Solver(mesh, opts){initParameters();}
void Solver::initParameters()
{
  _parameters.strings["application"] = "linear";
  _parameters.strings["method"] = "nitsche";
  _parameters.doubles["kin"] = 1.0;
  _parameters.doubles["kex"] = 1.0;
  _parameters.doubles["gamma"] = 2.0;
  _parameters.doubles["xgamma"] = 0.2;
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
std::unique_ptr<solvers::MeshUnitWithDataInterface> Solver::newMeshUnitWithDataPlain(const mesh::MeshUnitInterface* mesh) const
{
  return std::unique_ptr<solvers::MeshUnitWithDataInterface>(new MeshUnitWithData(_parameters.strings["application"]));
}
/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::ApplicationInterface> Solver::newApplication() const
{
  return std::unique_ptr<solvers::ApplicationInterface>(new Application(_parameters.strings["application"]));
}
std::unique_ptr<solvers::ModelInterface> Solver::newModel() const
{
  return std::unique_ptr<solvers::ModelInterface>(new Model());
}
/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::FemInterface> Solver::newFem(const std::string& varname, solverEnums::fem::femtype fem) const
{
  if(varname=="Lambda")
  {
    assert(_parameters.strings["method"] == "strong");
    return std::unique_ptr<solvers::FemInterface>(new P1Interface());
  }
  if(fem==solverEnums::fem::OWN)
  {
    if(_parameters.strings["method"] == "nitsche"
    or _parameters.strings["method"] == "strong")
    {
      return std::unique_ptr<solvers::FemInterface>(new P1Hansbo());
    }
    else if(_parameters.strings["method"] == "nitschestabc"
    or _parameters.strings["method"] == "nitschestabnc"
    or _parameters.strings["method"] == "newnitsche"
    or _parameters.strings["method"] == "newnitsche2")
    {
      return std::unique_ptr<solvers::FemInterface>(new P1Cut());
    }
    else if(_parameters.strings["method"] == "P1IFEM")
    {
      return std::unique_ptr<solvers::FemInterface>(new P1IIM());
    }
    else if(_parameters.strings["method"] == "CR1IFEM")
    {
      return std::unique_ptr<solvers::FemInterface>(new CR1IIM());
    }
  }
  return solvers::Solver::newFem(varname, fem);
}

/*--------------------------------------------------------------------------*/
void Solver::defineVariablesAndPdeParts()
{
  addVariablePlain("U", 1, solverEnums::fem::OWN);
  addDataVariablePlain("Phi", 1, solverEnums::fem::P1);
  if(_parameters.strings["method"] == "strong")
  {
    addVariablePlain("Lambda", 1, solverEnums::fem::OWN);
    alat::StringList vars("U,Lambda");
    addPdePartPlain("multiplier", std::unique_ptr<solvers::PdePartInterface>(new PdePartMultiplier(vars)));
  }
  alat::StringList vars("U");
  if(_parameters.strings["method"] == "P1IFEM" or _parameters.strings["method"] == "CR1IFEM")
  {
    addPdePartPlain("diffusion", std::unique_ptr<solvers::PdePartInterface>(new PdePartIIM(vars)));
  }
  else
  {
    addPdePartPlain("diffusion", std::unique_ptr<solvers::PdePartInterface>(new PdePart(vars)));
    // addPdePartPlain("diffusion", std::unique_ptr<solvers::PdePartInterface>(new PdePartCutInterface(vars)));
  }
}
/*--------------------------------------------------------------------------*/
void Solver::init()
{
  std::cerr << "Solver::init()\n";
  solvers::Solver::init();
  setOutputOptions(solver_options::solutionfilename) = _parameters.strings["method"];
}
/*--------------------------------------------------------------------------*/
solvers::ErrorsMap Solver::run()
{
  std::cerr << "Solver::run()\n";
  std::shared_ptr<MeshUnitWithData> meshunitwithdata = std::dynamic_pointer_cast<MeshUnitWithData>(_plainmeshunitwithdata);
  assert(meshunitwithdata);

  saveData();

  solvers::Solver::run();


  if(_parameters.strings["method"] == "nitsche"
  or _parameters.strings["method"] == "nitschestabc"
  or _parameters.strings["method"] == "nitschestabnc"
  or _parameters.strings["method"] == "strong"
  or _parameters.strings["method"] == "newnitsche"
  or _parameters.strings["method"] == "newnitsche2")
  {
    std::string filename = _output_manager.getSolutionFileName();
    meshunitwithdata->writeVtkCut("U", filename);
  }
  else if(_parameters.strings["method"] == "P1IFEM")
  {
    std::string filename = _output_manager.getSolutionFileName();
    meshunitwithdata->writeVtkIIM("U", filename);
  }

  bool test = _plainmeshunitwithdata->getApplication()->hasExactSolution("U");
  std::cerr << "test="<<test<<"\n";
  if(test)
  {
    solver_options::errors::opts opts = solver_options::errors::H1;
    opts += solver_options::errors::L2;
    // opts += solver_options::errors::L1;
    opts += solver_options::errors::Linf;
    opts += solver_options::errors::E;
    // std::cerr << "ops: " << opts << "\n";
    solvers::ErrorsMap errors = computeErrors("U", opts);
    // std::cerr << "errors: " << errors << "\n";
    return errors;
  }
  else
  {
    return solvers::ErrorsMap();
  }
}

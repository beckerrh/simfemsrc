#include  "Alat/directoryandfiles.hpp"
#include  "Alat/ghostmatrix.hpp"
#include  "Alat/ghostvector.hpp"
#include  "Alat/matrixonevariablearma.hpp"
#include  "Alat/matrixonevariable.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Mesh/boundaryinformation.hpp"
#include  "Mesh/mesh.hpp"
#include  "Perulangan/newtonsimple.hpp"
#include  "Perulangan/newtonlinesearch.hpp"
#include  "Solvers/application.hpp"
#include  "Solvers/p1.hpp"
#include  "Solvers/p2.hpp"
#include  "Solvers/cr1.hpp"
#include  "Solvers/rt0.hpp"
#include  "Solvers/solver.hpp"
#include  "Solvers/xdmfwriter.hpp"
#include  "Solvers/meshunitwithdata.hpp"
#include  "Solvers/nonlinearsolvervisitor.hpp"
#include  "Solvers/pdepartp1.hpp"
#include  <cassert>
#include  <sstream>
#include  <iomanip>

using namespace solvers;

/*--------------------------------------------------------------------------*/
Solver::~Solver()
{
  if(_printlevel) _chronometer.print(std::cerr);
}
Solver::Solver( const Solver& solver) : solvers::SolverInterface(solver), _opts(solver_options::none)
{
  assert(0);
}
Solver::Solver(std::shared_ptr<mesh::MeshInterface> mesh, solver_options::opts opts, int printlevel) : solvers::SolverInterface(), _mesh(mesh), _opts(opts), _partion_id(1), _initcalled(false), _meshsaved(false), _printlevel(printlevel)
{
  if(mesh) {_partion_id = mesh->getPartionId();}
  _parameters.strings["fem"] = "P1";
}

Solver& Solver::operator=( const Solver& solver)
{
  assert(0);
  solvers::SolverInterface::operator=(solver);
  return *this;
}
std::string Solver::getClassName() const
{
  return "Solver";
}
Solver* Solver::clone() const
{
  return new Solver(*this);
}
std::string& Solver::setOutputOptions(solver_options::output_manager_data opt)
{
  return _output_manager.data[opt];
}
void Solver::setTimeMesh(const mesh::TimeMeshData timemeshdata)
{
  _timemesh.setData(timemeshdata);
}
void  Solver::setSolverData(const solvers::solverdata solverdata)
{
  _solverdata = solverdata;
}
perulangan::NonlinearSolverInterface* Solver::getNonlinearSolver() {return _nonlinearsolver.get();}
const perulangan::NonlinearSolverInterface* Solver::getNonlinearSolver() const {return _nonlinearsolver.get();}

void Solver::setMesh(std::shared_ptr<mesh::MeshInterface> mesh) {_mesh=mesh;}
void Solver::setParameter(std::string name, std::string value)
{
  if(not _parameters.strings.hasKey(name)) {_error_string("setParameter","no such parameter",name);}
  _parameters.strings[name] = value;
}
void Solver::setParameter(std::string name, int value)
{
  if(not _parameters.ints.hasKey(name)) {_error_string("setParameter","no such parameter",name);}
  _parameters.ints[name] = value;
}
void Solver::setParameter(std::string name, double value)
{
  if(not _parameters.doubles.hasKey(name)) {_error_string("setParameter","no such parameter",name);}
  _parameters.doubles[name] = value;
}
void Solver::setParameter(std::string name, bool value)
{
  if(not _parameters.bools.hasKey(name)) {_error_string("setParameter","no such parameter",name);}
  _parameters.bools[name] = value;
}
/*---------------------------------------------------------*/
std::string Solver::getMeshInfo() const
{
  return _mesh->getInfo();
}
/*---------------------------------------------------------*/
std::string Solver::getInfo() const
{
  assert(_initcalled);
  std::string sep("-----------------\n");
  std::stringstream ss;
  ss << sep << getClassName() << " opts=" << _opts << "\n";
  ss << "parameters (sring): " << _parameters.strings << "\n";
  ss << "parameters (double): " << _parameters.doubles << "\n";
  ss << "parameters (int): " << _parameters.ints << "\n";
  ss << "parameters (bool): " << _parameters.bools << "\n";
  ss << _output_manager;
  ss << "Plain " <<  _plainmeshunitwithdata->getInfo();
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin();p!=_boundarymeshunitswithdata.end();p++)
  {
    ss << "Boundary" << p->first << " : " << p->second->getInfo();
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin();p!=_interfacemeshunitswithdata.end();p++)
  {
    ss << "Interface" << p->first << " : " << p->second->getInfo();
  }
  ss << sep;
  return ss.str();
}
/*--------------------------------------------------------------------------*/
std::ostream& solvers::operator<<(std::ostream &os, const Solver& solver)
{
  os << solver.getInfo() << "\n";
  return os;
}

/*---------------------------------------------------------*/
void Solver::writeXdmf() const
{
  std::string meshfilename  = _output_manager.getMeshFileName();
  if(not _meshsaved)
  {
    _mesh->saveH5(meshfilename+".h5");
    _meshsaved=true;
  }
  std::string filenamexdmf = _output_manager.getSolutionFileName() + ".xdmf";
  XdmfWriter xdmfwriter(&_output_manager, &_timemesh, _mesh->getPlainMesh());
  xdmfwriter.writeSolution(filenamexdmf, "u", _plainmeshunitwithdata->getVars());
}
/*---------------------------------------------------------*/
void Solver::writeDataXdmf(int it) const
{
  std::string meshfilename  = _output_manager.getMeshFileName();
  if(not _meshsaved)
  {
    _mesh->saveH5(meshfilename+".h5");
    _meshsaved=true;
  }
  std::string filenamexdmf = _output_manager.getDataFileName() + ".xdmf";
  // std::cerr << "filenamexdmf="<<filenamexdmf<<"\n";
  // assert(0);
  XdmfWriter xdmfwriter(&_output_manager, &_timemesh, _mesh->getPlainMesh());
  xdmfwriter.writeData(filenamexdmf, _plainmeshunitwithdata->getVarsData(), it);
}

/*--------------------------------------------------------------------------*/
void Solver::loadMesh(std::string meshtype, std::string meshname)
{
  // std::cerr << "Solver::loadMesh() meshtype="<<meshtype<<"\n";
  // assert(_mesh==nullptr);
  // setOutputOptions(solver_options::meshfilename) = meshname;
  // std::string meshfilename  = _output_manager.getMeshFileName();
  _mesh = mesh::Mesh::create(meshtype, _partion_id);
  // _mesh->setPartionId(_partion_id);
  // std::cerr << "Solver::loadMesh() meshfilename="<<meshfilename<<"\n";
  _mesh->loadH5(meshname);
}
/*--------------------------------------------------------------------------*/
void Solver::loadSolution(std::string filename)
{
  assert(0);
  _plainmeshunitwithdata->loadSolution(alat::GhostVector("u"), filename);
}
/*--------------------------------------------------------------------------*/
void Solver::saveSolution(const alat::GhostVector ghost, int it) const
{
  _chronometer.start("saveSolution");
  std::string filename = _output_manager.getSolutionFileName(meshEnums::Plain, ghost.getClassName(), it);
  _plainmeshunitwithdata->saveSolution(ghost, filename);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    std::stringstream ss;
    ss << _output_manager.getSolutionFileName(meshEnums::Boundary, ghost.getClassName(), it)<< p->first;
    p->second->saveSolution(ghost, ss.str());
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    std::stringstream ss;
    ss << _output_manager.getSolutionFileName(meshEnums::Interface, ghost.getClassName(), it)<< p->first;
    p->second->saveSolution(ghost, ss.str());
  }
  _chronometer.stop("saveSolution");
}
/*--------------------------------------------------------------------------*/
void Solver::saveData(int it) const
{
  _chronometer.start("saveSolution");
  std::string filename = _output_manager.getDataFileName(meshEnums::Plain, it);
  _plainmeshunitwithdata->saveData(filename);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    std::stringstream ss;
    ss << _output_manager.getDataFileName(meshEnums::Boundary, it)<< p->first;
    p->second->saveData(ss.str());
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    std::stringstream ss;
    ss << _output_manager.getDataFileName(meshEnums::Interface, it)<< p->first;
    p->second->saveData(ss.str());
  }
  _chronometer.stop("saveSolution");
}
/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::FemInterface> Solver::newFem(const std::string& varname, solverEnums::fem::femtype fem) const
{
  if(fem == solverEnums::fem::P1)
  {
    return std::unique_ptr<solvers::FemInterface>(new solvers::P1());
  }
  else if(fem == solverEnums::fem::RT0)
  {
    return std::unique_ptr<solvers::FemInterface>(new solvers::RT0());
  }
  else if(fem == solverEnums::fem::CR1)
  {
    return std::unique_ptr<solvers::FemInterface>(new solvers::CR1());
  }
  else if(fem == solverEnums::fem::P2)
  {
    return std::unique_ptr<solvers::FemInterface>(new solvers::P2());
  }
  _error_string("newFem", "unknwon fem", solverEnums::fem::femtypeToString(fem));
  return std::unique_ptr<solvers::FemInterface>(nullptr);
}
std::unique_ptr<solvers::ApplicationInterface> Solver::newApplication() const
{
  return std::unique_ptr<solvers::ApplicationInterface>(new solvers::Application());
}
std::unique_ptr<solvers::ModelInterface> Solver::newModel() const
{
  return std::unique_ptr<solvers::ModelInterface>(nullptr);
}

/*--------------------------------------------------------------------------*/
void Solver::defineVariablesAndPdeParts()
{
  addVariablePlain("U", 1, solverEnums::fem::P1);

  std::unique_ptr<solvers::PdePartInterface> pdepart;
  pdepart = std::unique_ptr<solvers::PdePartInterface>(new solvers::PdePartP1("U"));
  addPdePartPlain("laplace", std::move(pdepart));
  // std::unique_ptr<solvers::ApplicationInterface> application;
  // application = std::unique_ptr<solvers::ApplicationInterface>(new solvers::Application());
  // _plainmeshunitwithdata->setApplication(std::move(application));
}
// std::unique_ptr<solvers::PdePartInterface> Solver::newPdePart(const solvers::Variable& var) const
// {
//   std::unique_ptr<solvers::PdePartInterface> pdepart;
//   pdepart = std::unique_ptr<solvers::PdePartInterface>(new solvers::PdePartP1());
//   std::string fem = _parameters.strings["fem"];
//   if(fem=="P1")
//   {
//     pdepart->setFem(std::unique_ptr<solvers::FemInterface>(new solvers::P1()));
//   }
//   else
//   {
//     _error_string("newPdePart","unknown fem", fem);
//   }
//   return std::move(pdepart);
// }
/*--------------------------------------------------------------------------*/
std::unique_ptr<alat::MatrixOneVariableInterface> Solver::newMatrix() const
{
  if(solver_options::has_armamat(_opts))
  {
    return std::unique_ptr<alat::MatrixOneVariableInterface>(new alat::MatrixOneVariableArma);
  }
  else
  {
    return std::unique_ptr<alat::MatrixOneVariableInterface>(new alat::MatrixOneVariable);
  }
}
std::unique_ptr<alat::VectorOneVariableInterface> Solver::newVector() const
{
  return std::unique_ptr<alat::VectorOneVariableInterface>(new alat::VectorOneVariable);
}
std::unique_ptr<solvers::MeshUnitWithDataInterface> Solver::newMeshUnitWithDataPlain(const mesh::MeshUnitInterface* mesh) const
{
  return std::unique_ptr<solvers::MeshUnitWithDataInterface>(new MeshUnitWithData());
}
std::unique_ptr<solvers::MeshUnitWithDataInterface> Solver::newMeshUnitWithDataBoundary(int color, const mesh::MeshUnitInterface* mesh) const
{
  return std::unique_ptr<solvers::MeshUnitWithDataInterface>();
  // return std::unique_ptr<solvers::MeshUnitWithDataInterface>(new MeshUnitWithData);
}
std::unique_ptr<solvers::MeshUnitWithDataInterface> Solver::newMeshUnitWithDataInterface(int color, const mesh::MeshUnitInterface* mesh) const
{
  return std::unique_ptr<solvers::MeshUnitWithDataInterface>();
  // return std::unique_ptr<solvers::MeshUnitWithDataInterface>(new MeshUnitWithData);
}
std::unique_ptr<perulangan::NonlinearSolverInterface> Solver::newNonlinearSolver() const
{
  if(_solverdata.newton=="newton")
  {
    return std::unique_ptr<perulangan::NonlinearSolverInterface>(new perulangan::NewtonSimple);
  }
  else if(_solverdata.newton=="newtonls")
  {
    return std::unique_ptr<perulangan::NonlinearSolverInterface>(new perulangan::NewtonLineSearch("monotonicty"));
  }
  else if(_solverdata.newton=="newtonarmijo")
  {
    return std::unique_ptr<perulangan::NonlinearSolverInterface>(new perulangan::NewtonLineSearch("armijo"));
  }
  else if(_solverdata.newton=="newtonwolfe")
  {
    return std::unique_ptr<perulangan::NonlinearSolverInterface>(new perulangan::NewtonLineSearch("wolfe"));
  }
  else
  {
    _error_string("newNonlinearSolver","unknown nonlinear",_solverdata.newton);
    return std::unique_ptr<perulangan::NonlinearSolverInterface>(nullptr);
  }
}

/*--------------------------------------------------------------------------*/
void Solver::addDataVariablePlain(std::string name, int ncomp, solverEnums::fem::femtype fem)
{
  if(not _mesh) _error_string("addDataVariablePlain","no mesh");
  // _varinfodataplain.insert(varinfo(name, ncomp, fem));
  _varinfodataplain.push_back(varinfo(name, ncomp, fem));
}

/*--------------------------------------------------------------------------*/
void Solver::addVariablePlain(std::string name, int ncomp, solverEnums::fem::femtype fem)
{
  if(not _mesh) _error_string("addVariablePlain","no mesh");
  // _varinfoplain.insert(varinfo(name, ncomp, fem));
  _varinfoplain.push_back(varinfo(name, ncomp, fem));
  _var2ncomp[name] = ncomp;
}
/*--------------------------------------------------------------------------*/
void Solver::addPdePartPlain(std::string name, std::shared_ptr<solvers::PdePartInterface> pdepart)
{
  _plainmeshunitwithdata->addPdePart(name, pdepart);
}

/*--------------------------------------------------------------------------*/
void Solver::addBoundaryUnits(int first, int last)
{
  if(_initcalled) _error_string("init","init called already");
  if(not _mesh) _error_string("init","no mesh");
  if(first==std::numeric_limits<int>::max() and last==std::numeric_limits<int>::lowest())
  {
    const mesh::MeshInterface::BoundaryMeshUnitsMap& boundarymeshunitsmap = _mesh->getBoundaryMeshUnitsMap();
    const alat::IntSet keys = boundarymeshunitsmap.keys();
    _bdryunitscolors.insert(keys.begin(), keys.end());
  }
  else
  {
    for(int color=first; color<last;color++) _bdryunitscolors.insert(color);
  }
}
/*--------------------------------------------------------------------------*/
void Solver::enrolVector(alat::GhostVector ghost)
{
  _plainmeshunitwithdata->enrolVector(ghost);
}
void Solver::enrolMatrix(alat::GhostMatrix ghost)
{
  _plainmeshunitwithdata->enrolMatrix(ghost);
}
void Solver::initDataDir()
{
  std::string datadir = _output_manager.data[solver_options::datadir];
  std::string datadirold = datadir + ".old";
  // _path =   alat::_getPath() + "/";
  // std::cerr << "_path="<<_path<<"\n";
  std::string command;
  if(solver_options::has_erasedatadir(_opts))
  {
    if(alat::_directoryExists(datadirold))
    {
      command = "rm -rf " + datadirold;
      system( command.c_str() );
    }
    if(alat::_directoryExists(datadir))
    {
      command = "mv -f "+datadir + " " + datadirold;
      system( command.c_str() );
    }
  }
  if(not alat::_directoryExists(datadir))
  {
    command = "mkdir -p "+datadir;
    system( command.c_str() );
  }
  if(not alat::_directoryExists(datadir + "/Plain"))
  {
    command = "mkdir -p "+datadir + "/Plain";
    system( command.c_str() );
  }
  if(not alat::_directoryExists(datadir + "/Boundary"))
  {
    command = "mkdir -p "+datadir + "/Boundary";
    system( command.c_str() );
  }
  if(not alat::_directoryExists(datadir+ "/Interface"))
  {
    command = "mkdir -p "+datadir + "/Interface";
    system( command.c_str() );
  }
}
/*--------------------------------------------------------------------------*/
void Solver::defineGhosts()
{
  enrolVector(alat::GhostVector("u"));
  enrolVector(alat::GhostVector("f"));
  alat::GhostMatrix ghostmatrix("A");
  enrolMatrix(ghostmatrix);
  _plainmeshunitwithdata->enrolLinearSolver(alat::GhostLinearSolver("B","direct",ghostmatrix));
}

/*--------------------------------------------------------------------------*/
void Solver::init()
{
  // if(_initcalled) _error_string("init","init called already");

  if(_initcalled)
  {
    // std::cerr << "#########################\n";
    _chronometer.reset();
    _varinfoplain.clear(); _varinfodataplain.clear();
    _varinfobdry.clear(); _varinfoitfc.clear();
    _var2ncomp.clear();
    _bdryunitscolors.clear();
  }
  else
  {
    _chronometer.enrol("computeJacobian");
    _chronometer.enrol("computeLinearSolver");
    _chronometer.enrol("residual");
    _chronometer.enrol("solveLinear");
    _chronometer.enrol("solveNonlinear", false);
    _chronometer.enrol("saveSolution");
    _chronometer.enrol("errors");
  }

  // solver
  _nonlinearsolver = newNonlinearSolver();
  _nonlinearsolver->setVisitorPointer(
    std::unique_ptr<perulangan::NonlinearSolverVisitorInterface>(new NonlinearSolverVisitor(this))
  );
  // std::cerr << "Solver::init()  visitor = " << _nonlinearsolver->getVisitor()->getClassName() << "\n";

  // mesh
  assert(_mesh);
  _output_manager.init();
  initDataDir();
  const mesh::MeshInterface::BoundaryMeshUnitsMap& boundarymeshunitsmap = _mesh->getBoundaryMeshUnitsMap();
  const mesh::InterfaceMeshUnitsMap& interfacemeshunitsmap = _mesh->getInterfaceMeshUnitsMap();
  _plainmeshunitwithdata = newMeshUnitWithDataPlain(_mesh->getPlainMesh());
  for(alat::IntSet::const_iterator p=_bdryunitscolors.begin();p!=_bdryunitscolors.end();p++)
  {
    assert(boundarymeshunitsmap.hasKey(*p));
    _boundarymeshunitswithdata[*p] = newMeshUnitWithDataBoundary(*p, boundarymeshunitsmap[*p].get());
  }
  for(mesh::InterfaceMeshUnitsMap::const_iterator p=interfacemeshunitsmap.begin(); p!=interfacemeshunitsmap.end();p++)
  {
    _interfacemeshunitswithdata[p->first] = newMeshUnitWithDataInterface(p->first, p->second.get());
  }

  defineVariablesAndPdeParts();
  defineGhosts();

  assert(_plainmeshunitwithdata);
  // for(alat::Set<varinfo>::const_iterator p=_varinfoplain.begin(); p!=_varinfoplain.end();p++)
  for(alat::List<varinfo>::const_iterator p=_varinfoplain.begin(); p!=_varinfoplain.end();p++)
  {
    _plainmeshunitwithdata->addVariable(p->name, p->ncomp, p->femname);
  }
  for(alat::Map<int, varinfo>::const_iterator p=_varinfobdry.begin(); p!=_varinfobdry.end();p++)
  {
    _boundarymeshunitswithdata[p->first]->addVariable(p->second.name, p->second.ncomp, p->second.femname);
  }
  for(alat::Map<int, varinfo>::const_iterator p=_varinfoitfc.begin(); p!=_varinfoitfc.end();p++)
  {
    _interfacemeshunitswithdata[p->first]->addVariable(p->second.name, p->second.ncomp, p->second.femname);
  }
  // for(alat::Set<varinfo>::const_iterator p=_varinfodataplain.begin(); p!=_varinfodataplain.end();p++)
  for(alat::List<varinfo>::const_iterator p=_varinfodataplain.begin(); p!=_varinfodataplain.end();p++)
  {
    _plainmeshunitwithdata->addDataVariable(p->name, p->ncomp, p->femname);
  }

  _nonlinearsolver->init();

  _plainmeshunitwithdata->initMeshAndApplication(_mesh->getPlainMesh(), this);
  _plainmeshunitwithdata->initFemAndMemoryAndDataAndPdeParts();

  // for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  // {
  //   p->second->init(boundarymeshunitsmap[p->first].get(), this);
  // }
  // for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  // {
  //   p->second->init(interfacemeshunitsmap[p->first].get(), this);
  // }
  initData();
  saveData();
  writeDataXdmf(0);
  _initcalled = true;
}
/*--------------------------------------------------------------------------*/
void Solver::residual(perulanganEnums::residualstatus& status, alat::GhostVector& gr, const alat::GhostVector& gu, const alat::GhostVector& gf) const
{
  _chronometer.start("residual");
  status = perulanganEnums::ResidualStatusOk;
  _plainmeshunitwithdata->residual(status, gr, gu, gf);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->residual(status, gr, gu, gf);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->residual(status, gr, gu, gf);
  }
  _chronometer.stop("residual");
}
/*--------------------------------------------------------------------------*/
void Solver::initSolution(alat::GhostVector gu) const
{
  _plainmeshunitwithdata->initSolution(gu);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->initSolution(gu);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->initSolution(gu);
  }
}
/*--------------------------------------------------------------------------*/
void Solver::initData() const
{
  _plainmeshunitwithdata->initData();
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->initData();
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->initData();
  }
}
/*--------------------------------------------------------------------------*/
void Solver::setTimeMeshData(double time, double dt)
{
  _plainmeshunitwithdata->setTimeMeshData(time, dt);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->setTimeMeshData(time, dt);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->setTimeMeshData(time, dt);
  }
}
/*--------------------------------------------------------------------------*/
solvers::ErrorsMap Solver::computeErrors(std::string varname, solver_options::errors::opts opts) const
{
  _chronometer.start("errors");
  alat::GhostVector gu("u");
  solvers::ErrorsMap errorsmap;
  int ncomp = _var2ncomp[varname];
  if(solver_options::errors::has_L2(opts)) {errorsmap["L2"].set_size(ncomp);}
  if(solver_options::errors::has_H1(opts)) {errorsmap["H1"].set_size(ncomp);}
  if(solver_options::errors::has_L1(opts)) {errorsmap["L1"].set_size(ncomp);}
  if(solver_options::errors::has_E(opts)) {errorsmap["E"].set_size(ncomp);}
  if(solver_options::errors::has_Linf(opts)) {errorsmap["Linf"].set_size(ncomp);}
  for(solvers::ErrorsMap::iterator p=errorsmap.begin();p!=errorsmap.end();p++) {p->second.zeros();}
  _plainmeshunitwithdata->computeErrors(errorsmap, varname, gu);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->computeErrors(errorsmap, varname, gu);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->computeErrors(errorsmap, varname, gu);
  }
  if(solver_options::errors::has_L2(opts)) {errorsmap["L2"]=sqrt(errorsmap["L2"]);}
  if(solver_options::errors::has_H1(opts)) {errorsmap["H1"]=sqrt(errorsmap["H1"]);}
  if(solver_options::errors::has_E(opts)) {errorsmap["E"]=sqrt(errorsmap["E"]);}
  _chronometer.stop("errors");
  return errorsmap;
}

/*--------------------------------------------------------------------------*/
void Solver::computeRhs(alat::GhostVector gf, const alat::GhostVector gu) const
{
  _plainmeshunitwithdata->computeRhs(gf, gu);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->computeRhs(gf, gu);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->computeRhs(gf, gu);
  }
}
/*--------------------------------------------------------------------------*/
void Solver::constructMatrixAndLinearSolvers(perulanganEnums::matrixstatus& status, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u)
{
  _chronometer.start("computeJacobian");
  _plainmeshunitwithdata->computeJacobian(status, A, u);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->computeJacobian(status, A, u);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->computeJacobian(status, A, u);
  }
  _chronometer.stop("computeJacobian");
  _chronometer.start("computeLinearSolver");
  _plainmeshunitwithdata->computeLinearSolver(status, B, A, u);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->computeLinearSolver(status, B, A, u);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->computeLinearSolver(status, B, A, u);
  }
  _chronometer.stop("computeLinearSolver");
}
/*--------------------------------------------------------------------------*/
int Solver::solveLinear(perulanganEnums::iterationstatus& status, const alat::GhostLinearSolver& B, const alat::GhostMatrix& A, alat::GhostVector& du, const alat::GhostVector& r) const
{
  _chronometer.start("solveLinear");
  int niter = _plainmeshunitwithdata->solveLinear(status, B, A, du, r);
  _chronometer.stop("solveLinear");
  return niter;
}
/*--------------------------------------------------------------------------*/
void Solver::setLinearTolerance(double rtol, double gtol, alat::GhostLinearSolver& B)
{
  _plainmeshunitwithdata->getLinearSolver(B)->setTolerance(rtol, gtol);
}
/*--------------------------------------------------------------------------*/
void Solver::vectorZero(alat::GhostVector& gu) const
{
  _plainmeshunitwithdata->vectorZero(gu);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->vectorZero(gu);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->vectorZero(gu);
  }
}
void Solver::vectorAdd(alat::GhostVector& gu, double s, const alat::GhostVector& gv) const
{
  _plainmeshunitwithdata->vectorAdd(gu, s, gv);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->vectorAdd(gu, s, gv);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->vectorAdd(gu, s, gv);
  }
}
double Solver::vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const
{
  double dot = _plainmeshunitwithdata->vectorDot(gu, gv);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    dot += p->second->vectorDot(gu, gv);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    dot += p->second->vectorDot(gu, gv);
  }
  return dot;
}
void Solver::vectorEqual(alat::GhostVector& gu, const alat::GhostVector& gv) const
{
  _plainmeshunitwithdata->vectorEqual(gu, gv);
  for(MeshUnitWithDataMap::const_iterator p=_boundarymeshunitswithdata.begin(); p!=_boundarymeshunitswithdata.end();p++)
  {
    p->second->vectorEqual(gu, gv);
  }
  for(MeshUnitWithDataMap::const_iterator p=_interfacemeshunitswithdata.begin(); p!=_interfacemeshunitswithdata.end();p++)
  {
    p->second->vectorEqual(gu, gv);
  }
}

/*--------------------------------------------------------------------------*/
void Solver::run()
{
  perulangan::NewtonOutputData newtonoutputdata;
  perulanganEnums::residualstatus status;
  alat::GhostMatrix gA("A");
  alat::GhostVector gu("u");
  alat::GhostVector gf("f");
  alat::GhostLinearSolver gB("B","direct",gA);
  assert(_initcalled);
  initSolution(gu);
  double time = _timemesh.t(0);
  saveSolution(gu, 0);
  double  timestepratio = _solverdata.timestepratio;
  for(int it=1;it<_timemesh.n();it++)
  {
    double time_old = time;
    time = _timemesh.t(it);
    double dt = time-time_old;
    // if(it==1)
    // {
    //   setTimeMeshData(time, _solverdata.timestepratio*dt);
    //   constructMatrixAndLinearSolvers(status, gB, gA, gu);
    // }
    if(solver_options::has_dynamic(_opts))
    {
      int nmicrosteps = macroStep(time_old, time, timestepratio*dt, gA, gB, gu, gf, newtonoutputdata);
      std::cerr << "Solver::macroStep() " << it << " time="<<time << " dt="<<dt<<" nmicrosteps="<<nmicrosteps<<"\n";
      timestepratio = 1.0/nmicrosteps;
      // std::cerr << "Solver::macroStep() ?? " << 1.0/nmicrosteps << "\n";
    }
    else
    {
      computeRhs(gf, gu);
      _chronometer.start("solveNonlinear");
      _nonlinearsolver->solve(newtonoutputdata, gB, gA, gu, gf);
      // std::cerr << "netwondata="<<newtondata<<"\n";
      _chronometer.stop("solveNonlinear");
    }
    saveSolution(gu, it);
  }
  std::cerr << "ntotal new matrix " << newtonoutputdata.nredo_matrix << "\n";

}

/*--------------------------------------------------------------------------*/
int Solver::macroStep(double timebegin, double timeend, double dtfirst, alat::GhostMatrix gA, alat::GhostLinearSolver gB, alat::GhostVector gu, alat::GhostVector gf, perulangan::NewtonOutputData& newtonoutputdata)
{
  assert(timebegin+dtfirst<=timeend);
  double time = timebegin;
  double dt = dtfirst;
  int nstepsnomat=0;
  int nsteps=0;
  bool dynamic=true;
  if(dynamic)
  {
    perulangan::NewtonInputData newtoninputdata(_nonlinearsolver->getNewtonInputData());
    newtoninputdata.printstep=0;
    newtoninputdata.printlinesearch=false;
    getNonlinearSolver()->setNewtonInputData(newtoninputdata);
  }
  while(time<timeend)
  {
    nsteps++;
    time += dt;
    double tmissing = timeend-time;
    setTimeMeshData(time, dt);
    computeRhs(gf, gu);
    _chronometer.start("solveNonlinear");
    int nmat0 = newtonoutputdata.nredo_matrix;
    _nonlinearsolver->solve(newtonoutputdata, gB, gA, gu, gf);
    _chronometer.stop("solveNonlinear");
    if(newtonoutputdata.newton_status != perulanganEnums::NewtonStatusConverged)
    {
      std::cerr << "newtonoutputdata="<<newtonoutputdata<<"\n";
      // redo timestep
      assert(0);
    }
    int nmat = newtonoutputdata.nredo_matrix-nmat0;
    std::stringstream info;
    info << " "<< std::setiosflags(std::ios::scientific) << std::setprecision(3) << std::setw(4)<< nsteps;
    info << " == " << time <<  "  ~"  << tmissing <<  "~ " << dt << " == ";
    info << newtonoutputdata.niter << "  *" << nmat << " " << nstepsnomat << "* " << "\n";
    std::cerr << info.str();
    if(time>=timeend) return nsteps;
    if(nmat>1)
    {
      nstepsnomat=0;
      dt *= 0.5;
      setTimeMeshData(time, dt);
      _nonlinearsolver->getVisitor()->constructMatrixAndLinearSolvers(newtonoutputdata.matrix_status, gB, gA, gu);
      if(dt<1e-10)
      {
        exit(100);
      }
    }
    else if(nmat==0)
    {
      if(nstepsnomat>3)
      {
        dt *= 1.1;
      }
      nstepsnomat++;
    }
    else
    {
      nstepsnomat=0;
    }
    int nleft = (int) (tmissing/dt);
    if(nleft==0)
    {
      dt = tmissing;
    }
    else
    {
      dt = tmissing/nleft;
    }
  }
  _error_string("macroStep","what happens ?");
  return -1;
  // std::cerr << "Solver::solve() " << it << " time="<<time << " dt="<<dt<<"\n";
  // if(it==1) constructMatrixAndLinearSolvers(status, gB, gA, gu);
}

#include  "Alat/matrixonevariablearma.hpp"
#include  "Alat/matrixonevariable.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Mesh/boundaryinformation.hpp"
#include  "Mesh/measureofcell.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Mesh/sidesandcells.hpp"
#include  "Mesh/normals.hpp"
#include  "Perulangan/directlinearsolver.hpp"
#include  "Solvers/p1.hpp"
#include  "Solvers/meshunitwithdata.hpp"
#include  "Solvers/xdmfwriter.hpp"
#include  <cassert>
#include  <sstream>

using namespace solvers;

/*--------------------------------------------------------------------------*/
MeshUnitWithData::~MeshUnitWithData() {}
MeshUnitWithData::MeshUnitWithData(): MeshUnitWithDataInterface(), _initcalled(false), _mesh(NULL), _solver(NULL) {}
MeshUnitWithData::MeshUnitWithData( const MeshUnitWithData& meshunitwithdata): MeshUnitWithDataInterface(meshunitwithdata)
{
  assert(0);
}

MeshUnitWithData& MeshUnitWithData::operator=( const MeshUnitWithData& meshunitwithdata)
{
  assert(0);
  MeshUnitWithDataInterface::operator=(meshunitwithdata);
  return *this;
}
std::string MeshUnitWithData::getClassName() const
{
  return "MeshUnitWithData";
}
const PdePartsMap& MeshUnitWithData::getPdeParts() const {return _pdeparts;}
PdePartsMap& MeshUnitWithData::getPdeParts() {return _pdeparts;}
const alat::Vector<Variable>& MeshUnitWithData::getVars() const {return _vars;}
const alat::Vector<Variable>& MeshUnitWithData::getVarsData() const {return _varsdata;}
void MeshUnitWithData::setTimeMeshData(double time, double dt)
{
  for(solvers::PdePartsMap::iterator p=_pdeparts.begin();p!=_pdeparts.end();p++)
  {
    p->second->setTimeMeshData(time, dt);
  }
}

/*--------------------------------------------------------------------------*/
void MeshUnitWithData::enrolVector(alat::GhostVector ghost) {_vectoragency.enrol(ghost);}
void MeshUnitWithData::enrolMatrix(alat::GhostMatrix ghost) {_matrixagency.enrol(ghost);}
void MeshUnitWithData::enrolLinearSolver(alat::GhostLinearSolver ghost) {_linearsolveragency.enrol(ghost);}
alat::VectorAllVariables* MeshUnitWithData::getVector(alat::GhostVector ghost) const {return _vectoragency(ghost).get();}
alat::MatrixAllVariables* MeshUnitWithData::getMatrix(alat::GhostMatrix ghost) const {return _matrixagency(ghost).get();}
perulangan::LinearSolverInterface* MeshUnitWithData::getLinearSolver(alat::GhostLinearSolver ghost) const{return _linearsolveragency(ghost).get();}
const solvers::FemInterface* MeshUnitWithData::getFemData(std::string name) const
{
  return _femsdata[_var2indexdata[name]].get();
}
alat::VectorOneVariable* MeshUnitWithData::getDataVector(std::string name) const
{
  alat::GhostVector ghost("Data", "data");
  alat::VectorAllVariables& data = *getVector(ghost);
  alat::VectorOneVariable* vector = dynamic_cast<alat::VectorOneVariable*>(data.get(_var2indexdata[name]));
  assert(vector);
  return vector;
}

/*--------------------------------------------------------------------------*/
const solvers::ApplicationInterface* MeshUnitWithData::getApplication() const {assert(_application);return _application.get();}
const solvers::ModelInterface* MeshUnitWithData::getModel() const {return _model.get();}
solvers::FemMap* MeshUnitWithData::getFemMap() const {return &_fems;}
solvers::FemMap* MeshUnitWithData::getFemExMap() const {return &_femsex;}
const mesh::MeshUnitInterface* MeshUnitWithData::getMesh() const {return _mesh;}

/*--------------------------------------------------------------------------*/
void MeshUnitWithData::addDataVariable(std::string name, int ncomp, solverEnums::fem::femtype fem)
{
  if(_initcalled) _error_string("init","init called already");
  _varsdata.push_back(Variable(name, ncomp, fem));

}
void MeshUnitWithData::addVariable(std::string name, int ncomp, solverEnums::fem::femtype fem)
{
  if(_initcalled) _error_string("init","init called already");
  _vars.push_back(Variable(name, ncomp, fem));

}
void MeshUnitWithData::addPdePart(std::string name, std::shared_ptr<solvers::PdePartInterface> pdepart)
{
  _pdeparts[name] = std::move(pdepart);
  // _pdeparts.insert( std::make_pair(name,pdepart));
}

/*--------------------------------------------------------------------------*/
std::string MeshUnitWithData::getInfo() const
{
  std::stringstream ss;
  ss << "*System* ";
  for(int ivar=0; ivar<_vars.size(); ivar++)
  {
    ss << _vars[ivar].getName() << " " << _vars[ivar].getNcomp() << " ";
    ss << _fems[ivar]->getInfo() <<  "\t";
  }
  ss << "\n*Application* " << _application->getInfo();
  if(_model)
  {
    ss << "\n*Model* " << _model->getInfo();
  }
  ss << "\n*Data* ";
  for(int ivar=0; ivar<_varsdata.size(); ivar++)
  {
    ss << _varsdata[ivar].getName() << " " << _varsdata[ivar].getNcomp() << " ";
    ss << _femsdata[ivar]->getInfo() <<  "\t";
  }
  ss << "\n*PdeParts* ";
  for(solvers::PdePartsMap::const_iterator p=_pdeparts.begin();p!=_pdeparts.end();p++)
  {
    ss << p->first << "|" << p->second->getInfo() << " ";
  }
  ss << "\n*Vectors*";
  for(alat::AgencyVector::iterator p=_vectoragency.begin();p!=_vectoragency.end();p++)
  {
    ss << p->first << "(" << p->second->size() << ") ";
  }
  ss << "\n";
  return ss.str();
}

/*--------------------------------------------------------------------------*/
void MeshUnitWithData::loadSolution(alat::GhostVector ghost, std::string filename)
{
  alat::VectorAllVariables& u = *getVector(ghost);
  int nvars = _vars.size();
  assert(nvars==u.size());
  alat::VectorAllVariables uin(nvars);
  for(int ivar=0; ivar<nvars; ivar++)
  {
    int ncompi = _vars[ivar].getNcomp();
    int ni = _meshinfo->nnodes;
    uin[ivar] = _solver->newVector();
    uin[ivar]->set_size(ni*ncompi);
    uin[ivar]->setNcompAndN(ncompi, ni);
  }
  uin.loadhdf5(filename, _varnames);
  for(int ivar=0; ivar<nvars; ivar++)
  {
    solverEnums::fem::femtype fem = _vars[ivar].getFem();
    if(fem == solverEnums::fem::P1)
    {
      u[ivar]->equal(uin[ivar].get());
    }
    else
    {
      _fems[ivar]->fromP1(u[ivar].get(), uin[ivar].get());
    }
  }
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::saveSolution(const alat::GhostVector ghost, std::string filename) const
{
  alat::VectorAllVariables& u = *getVector(ghost);
  int nvars = _vars.size();
  assert(nvars==u.size());
  alat::VectorAllVariables uout(nvars);
  for(int ivar=0; ivar<nvars; ivar++)
  {
    solverEnums::fem::femtype fem = _vars[ivar].getFem();
    int ncompi = _vars[ivar].getNcomp();
    int ni = _meshinfo->nnodes;
    uout[ivar] = _solver->newVector();
    uout[ivar]->set_size(ni*ncompi);
    uout[ivar]->setNcompAndN(ncompi, ni);
    if(not _fems[ivar]->canInterpolateToP1()) continue;
    if(fem == solverEnums::fem::P1)
    {
      uout[ivar]->equal(u.get(ivar));
    }
    else
    {
      _fems[ivar]->toP1(uout[ivar].get(), u[ivar].get());
    }
  }
  uout.savehdf5(filename, _varnames);
  // const alat::VectorOneVariable* uv = dynamic_cast<const alat::VectorOneVariable*>(uout.get(0));
  // std::cerr << "MeshUnitWithData::saveSolution u=" << *uv << "\n";
  // getVector(ghost)->savehdf5(filename, _varnames);
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::saveData(std::string filename) const
{
  alat::GhostVector ghost("Data", "data");
  alat::VectorAllVariables& u = *getVector(ghost);
  int nvars = _varsdata.size();
  assert(nvars==u.size());
  alat::VectorAllVariables uout(nvars);
  for(int ivar=0; ivar<nvars; ivar++)
  {
    solverEnums::fem::femtype fem = _varsdata[ivar].getFem();
    int ncompi = _varsdata[ivar].getNcomp();
    int ni = _meshinfo->nnodes;
    uout[ivar] = _solver->newVector();
    uout[ivar]->set_size(ni*ncompi);
    uout[ivar]->setNcompAndN(ncompi, ni);
    // std::cerr << "data="<<_varnamesdata[ivar]<<"\n";
    // std::cerr << "fem="<<_femsdata[ivar]->getClassName()<<"\n";
    if(fem == solverEnums::fem::P1)
    {
      uout[ivar]->equal(u.get(ivar));
    }
    else
    {
      if(_femsdata[ivar]->canInterpolateToP1())
      _femsdata[ivar]->toP1(uout[ivar].get(), u[ivar].get());
    }
  }
  uout.savehdf5(filename, _varnamesdata);
  // const alat::VectorOneVariable* uv = dynamic_cast<const alat::VectorOneVariable*>(uout.get(0));
  // std::cerr << "MeshUnitWithData::saveSolution u=" << *uv << "\n";
  // getVector(ghost)->savehdf5(filename, _varnames);
}

/*--------------------------------------------------------------------------*/
void MeshUnitWithData::initMeshAndApplication(mesh::MeshUnitInterface* mesh, const solvers::SolverInterface* solver)
{
  if(_initcalled) _error_string("init","init called already");
  _initcalled = true;

  _mesh = mesh;
  _solver = solver;

  assert(_mesh->geometryObjectExists(meshEnums::MeasureOfCell));
  std::shared_ptr<const mesh::MeasureOfCell> _measureofcell = std::dynamic_pointer_cast<const mesh::MeasureOfCell>(_mesh->getGeometryObject(meshEnums::MeasureOfCell));
  assert(_measureofcell);

  assert(_mesh->geometryObjectExists(meshEnums::Normals));
  std::shared_ptr<const mesh::Normals> _normals = std::dynamic_pointer_cast<const mesh::Normals>(_mesh->getGeometryObject(meshEnums::Normals));
  assert(_normals);

  _meshinfo = std::unique_ptr<MeshInfo>(new MeshInfo(
    _mesh->getNodesAndNodesOfCells().getNodes(),
    _mesh->getNodesAndNodesOfCells().getNodesOfCells(),
    _mesh->getSidesAndCells().getSidesOfCells(),
    _mesh->getSidesAndCells().getCellsOfSides(),
    _mesh->getEdgesAndCells().getEdgesOfCells(),
    _mesh->getSidesAndCells().getNodesOfSides(),
    _mesh->getEdgesAndCells().getNodesOfEdges(),
    _mesh->getEdgesAndCells().getLocalnodesOfEdgesOfCells(),
    _measureofcell->getMeasureOfCell(),
    _normals->getNormals(),
    _normals->getSigma(),
    (double) _mesh->getDimension(),
    _mesh->getNNodes(),
    _mesh->getNCells(),
    _mesh->getNSides(),
    _mesh->getNEdges(),
    _mesh->getNNodesPerCell(),
    _mesh->getNEdgesPerCell(),
    _mesh->getNSidesPerCell(),
    _mesh->getNNodesPerSide(),
    _mesh->getBoundaryInformationMap()
  ));

  // std::cerr << "MeshUnitWithData::init() _ncells=" << _ncells << "\n";

  // varnames
  int nvars = _vars.size();
  _varnames.set_size(nvars);
  alat::armaivec ncomps(nvars);
  for(int ivar=0; ivar<nvars; ivar++)
  {
    // _vars[ivar].init(mesh);
    _varnames[ivar] = _vars[ivar].getName();
    _var2index[_vars[ivar].getName()] = ivar;
    ncomps[ivar] = _vars[ivar].getNcomp();
  }
  // std::cerr << "MeshUnitWithData::init() _varnames=" << _varnames << "\n";
  // std::cerr << "MeshUnitWithData::init() _var2index=" << _var2index << "\n";

  // Data
  int nvarsdata = _varsdata.size();
  _varnamesdata.set_size(nvarsdata);
  alat::armaivec ncompsdata(nvarsdata);
  for(int ivar=0; ivar<nvarsdata; ivar++)
  {
    // _vars[ivar].init(mesh);
    _varnamesdata[ivar] = _varsdata[ivar].getName();
    _var2indexdata[_varsdata[ivar].getName()] = ivar;
    ncompsdata[ivar] = _varsdata[ivar].getNcomp();
  }
  // Model and Application
  _model = _solver->newModel();
  if(_model) _model->initModel(mesh, _solver->getParameters());
  assert(not _application);
  _application = _solver->newApplication();
  if(_model) _application->setModel(_model.get());
  _application->initApplication(mesh, _varnames, _varnamesdata, ncomps, ncompsdata);
}

/*--------------------------------------------------------------------------*/
void MeshUnitWithData::initFemAndMemoryAndDataAndPdeParts()
{
  int nvars = _vars.size();
  int nvarsdata = _varsdata.size();
  _ncomps.set_size(nvars);
  _nlocals.set_size(nvars);
  // fem
  _fems.set_size(nvars);
  _femsex.set_size(nvars);
  for(int ivar=0; ivar<nvars; ivar++)
  {
    _fems[ivar] = _solver->newFem(_varnames[ivar], _vars[ivar].getFem());
    _femsex[ivar] = _solver->newFem(_varnames[ivar], _vars[ivar].getFem());
    // std::cerr << "MeshUnitWithData::init() fem="<<_fems[ivar]->getInfo()<<"\n";
    int ncomp = _vars[ivar].getNcomp();
    _fems[ivar]->initFem(ivar, _mesh, _meshinfo.get(), ncomp);
    _femsex[ivar]->initFem(ivar, _mesh, _meshinfo.get(), ncomp);
    _ncomps[ivar] = ncomp;
    _nlocals[ivar] = _fems[ivar]->getNPerCell();
  }
  // Vectors
  for(alat::AgencyVector::iterator p=_vectoragency.begin();p!=_vectoragency.end();p++)
  {
    alat::VectorAllVariables& vector = *p->second.get();
    vector.set_size(nvars);
    for(int ivar=0; ivar<nvars; ivar++)
    {
      int ncompi = _vars[ivar].getNcomp();
      int ni = _fems[ivar]->getN();
      // std::cerr << "vector="<<p->first << " ni="<<ni<< "\n";
      vector[ivar] = _solver->newVector();
      vector[ivar]->set_size(ni*ncompi);
      vector[ivar]->setNcompAndN(ncompi, ni);
    }
  }
  // Data
  _femsdata.set_size(nvarsdata);
  for(int ivar=0; ivar<nvarsdata; ivar++)
  {
    _femsdata[ivar] = _solver->newFem(_varnamesdata[ivar], _varsdata[ivar].getFem());
    // std::cerr << "MeshUnitWithData::init() fem data ="<< _varnamesdata[ivar] << " " <<_femsdata[ivar]->getInfo()<<"\n";
    int ncomp = _varsdata[ivar].getNcomp();
    _femsdata[ivar]->initFem(ivar, _mesh, _meshinfo.get(), ncomp);
  }
  alat::GhostVector ghost("Data", "data");
  _vectoragency.enrol(ghost);
  std::shared_ptr<alat::VectorAllVariables> datavector = _vectoragency[ghost];
  datavector->set_size(nvarsdata);
  for(int ivar=0;ivar<nvarsdata;ivar++)
  {
    int ncompi = _varsdata[ivar].getNcomp();
    int ni = _femsdata[ivar]->getN();
    (*datavector)[ivar] = _solver->newVector();
    (*datavector)[ivar]->set_size(ni*ncompi);
    (*datavector)[ivar]->setNcompAndN(ncompi, ni);
  }

  // pdeparts
  _pdepartdata.set_size(_ncomps, _nlocals);
  for(PdePartsMap::iterator p=_pdeparts.begin();p!=_pdeparts.end();p++)
  {
    assert(p->second);
    p->second->initPdePart(_meshinfo.get(), _var2index, this, _solver->getParameters());
  }
  // Matrix
  for(alat::AgencyMatrix::iterator p=_matrixagency.begin();p!=_matrixagency.end();p++)
  {
    alat::MatrixAllVariables& matrix = *p->second.get();
    matrix.set_size(nvars, nvars);
    bool needsconnectivity=true;
    for(int ivar=0; ivar<nvars; ivar++)
    {
      int ncompi = _vars[ivar].getNcomp();
      int ni = _fems[ivar]->getN();
      for(int jvar=0; jvar<nvars; jvar++)
      {
        int ncompj = _vars[jvar].getNcomp();
        int nj = _fems[jvar]->getN();
        matrix(ivar,jvar) = _solver->newMatrix();
        matrix(ivar,jvar)->set_size(ni*ncompi, nj*ncompj);
        if(not matrix(ivar,jvar)->needsConnectivity())
        {
          needsconnectivity=false;
        }
        else
        {
          assert(needsconnectivity);
        }
      }
    }
    if(not needsconnectivity) {continue;}

    alat::Matrix<alat::SparsityPatternSoft> sparsitypatternsoft(nvars,nvars);
    for(int ivar=0; ivar<nvars; ivar++)
    {
      int ncompi = _vars[ivar].getNcomp();
      int ni = _fems[ivar]->getN();
      for(int jvar=0; jvar<nvars; jvar++)
      {
        sparsitypatternsoft(ivar,jvar).set_size(ni*ncompi);
      }
    }
    for(PdePartsMap::const_iterator p=_pdeparts.begin(); p!=_pdeparts.end();p++)
    {
      const alat::armaivec& ivars = p->second->getIvars();
      for(int ii=0; ii<ivars.size(); ii++)
      {
        int ivar = ivars[ii];
        int ncompi = _vars[ivar].getNcomp();
        int ni = _fems[ivar]->getN();
        for(int jj=0; jj<ivars.size(); jj++)
        {
          int jvar = ivars[jj];
          int ncompj = _vars[jvar].getNcomp();
          int nj = _fems[jvar]->getN();
          if(p->second->loopCells())
          {
            for(int iK=0; iK<_meshinfo->ncells;iK++)
            {
              _fems[ivar]->setCell(iK);
              _fems[jvar]->setCell(iK);
              int nloci = _fems[ivar]->getNPerCell(iK);
              int nlocj = _fems[jvar]->getNPerCell(iK);
              alat::armaivec vec_i(ncompi*nloci);
              alat::armaivec vec_j(ncompj*nlocj);
              _fems[ivar]->setVectorIndices(iK, vec_i, ncompi);
              _fems[jvar]->setVectorIndices(iK, vec_j, ncompj);
              for(alat::armaimat::const_iterator p=vec_i.begin();p!=vec_i.end();p++)
              {
                for(alat::armaimat::const_iterator q=vec_j.begin();q!=vec_j.end();q++)
                {
                  sparsitypatternsoft(ivar,jvar)[*p].insert(*q);
                }
              }
            }
          }
          // if(p->second->loopBoundary())
          if(p->second->loopInteriorSides())
          {
            for(int iS=0; iS<_meshinfo->nsides;iS++)
            {
              int iKin = _meshinfo->cells_of_sides(0,iS);
              int iKex = _meshinfo->cells_of_sides(1,iS);
              if(iKex<0) continue;
              if(not p->second->interiorsidecoupling(iKin, iKex)) continue;
              _fems[ivar]->setCell(iKin);
              _fems[jvar]->setCell(iKin);
              _femsex[ivar]->setCell(iKex);
              _femsex[jvar]->setCell(iKex);
              int nloci_in = _fems[ivar]->getNPerCell(iKin);
              int nlocj_in = _fems[jvar]->getNPerCell(iKin);
              int nloci_ex = _femsex[ivar]->getNPerCell(iKex);
              int nlocj_ex = _femsex[jvar]->getNPerCell(iKex);
              alat::armaivec vec_i_in(ncompi*nloci_in);
              alat::armaivec vec_j_in(ncompj*nlocj_in);
              alat::armaivec vec_i_ex(ncompi*nloci_ex);
              alat::armaivec vec_j_ex(ncompj*nlocj_ex);
              _fems[ivar]->setVectorIndices(iKin, vec_i_in, ncompi);
              _fems[jvar]->setVectorIndices(iKin, vec_j_in, ncompj);
              _femsex[ivar]->setVectorIndices(iKex, vec_i_ex, ncompi);
              _femsex[jvar]->setVectorIndices(iKex, vec_j_ex, ncompj);
              for(alat::armaimat::const_iterator p=vec_i_in.begin();p!=vec_i_in.end();p++)
              {
                for(alat::armaimat::const_iterator q=vec_j_ex.begin();q!=vec_j_ex.end();q++)
                {
                  sparsitypatternsoft(ivar,jvar)[*p].insert(*q);
                }
              }
              for(alat::armaimat::const_iterator p=vec_i_ex.begin();p!=vec_i_ex.end();p++)
              {
                for(alat::armaimat::const_iterator q=vec_j_in.begin();q!=vec_j_in.end();q++)
                {
                  sparsitypatternsoft(ivar,jvar)[*p].insert(*q);
                }
              }
            }
          }
        }
      }
      p->second->additionCouplings(sparsitypatternsoft);
    }
    for(int ivar=0; ivar<nvars; ivar++)
    {
      for(int jvar=0; jvar<nvars; jvar++)
      {
        matrix(ivar,jvar)->initSparsityPattern(sparsitypatternsoft(ivar,jvar));
      }
    }
  }
  for(alat::AgencyLinearSolver::iterator p=_linearsolveragency.begin();p!=_linearsolveragency.end();p++)
  {
    const std::string& name = p->first.getName();
    const std::string& type = p->first.getType();
    const alat::GhostMatrix& ghostmatrix = p->first.getMatrix();
    // std::cerr << "ghostmatrix=" <<ghostmatrix<<"\n";
    // std::cerr << "type=" <<type<<"\n";
    // std::cerr << "name=" <<name<<"\n";
    const alat::MatrixAllVariables* matrix = _matrixagency[ghostmatrix].get();
    if(type=="direct")
    {
      p->second = std::unique_ptr<perulangan::LinearSolverInterface>
      (
        new perulangan::DirectLinearSolver(matrix, this)
      );
    }
    else {assert(0);}
    p->second->reInit();
  }
}

/*--------------------------------------------------------------------------*/
void MeshUnitWithData::initSolution(alat::GhostVector gu) const
{
  alat::VectorAllVariables& u = *getVector(gu);
  u.fillzeros();
  int nvars = _vars.size();
  assert(nvars==u.size());
  const alat::VectorOneVariable* uv = dynamic_cast<const alat::VectorOneVariable*>(u.get(0));
  for(int ivar=0; ivar<nvars; ivar++)
  {
    if(_application->hasInitialCondition(ivar))
    {
      _fems[ivar]->interpolate(u.get(ivar), _application->getInitialCondition(ivar));
      // std::cerr << "#1 MeshUnitWithData::initSolution u=" << *uv << "\n";
    }
    else
    {
      u.get(ivar)->equal(0.0);
      // std::cerr << "#1 zero MeshUnitWithData::initSolution u=" << *uv << "\n";
    }
    if(_application->hasDirichlet(ivar))
    {
      // std::cerr << "#2 IN MeshUnitWithData::initSolution u=" << *uv << "\n";
      _fems[ivar]->strongDirichlet(u.get(ivar), _application->getDirichlet(ivar), _application->getStrongDirichletColor());
      // std::cerr << "#2 OUT MeshUnitWithData::initSolution u=" << *uv << "\n";
    }
  }
}
/*-------------------------------------------------------------------------*/
void MeshUnitWithData::initDataVector(alat::VectorOneVariableInterface* u, int ivar, std::string varname) const
{
  if(not _application->hasDataFunction(ivar))
  {
    _error_string("initData", "application does not have DataFunction for", _varnamesdata[ivar]);
  }
  _femsdata[ivar]->interpolate(u, _application->getDataFunction(ivar));
}
/*-------------------------------------------------------------------------*/
void MeshUnitWithData::initData() const
{
  alat::GhostVector ghost("Data", "data");
  alat::VectorAllVariables& data = *getVector(ghost);
  data.fillzeros();
  int nvarsdata = _varsdata.size();
  assert(nvarsdata==data.size());
  for(int ivar=0; ivar<nvarsdata; ivar++)
  {
    initDataVector(data.get(ivar), ivar, _varnamesdata[ivar]);
    // if(not _application->hasDataFunction(ivar))
    // {
    //   _error_string("initData", "application does not have DataFunction for", _varnamesdata[ivar]);
    // }
    // _femsdata[ivar]->interpolate(data.get(ivar), _application->getDataFunction(ivar));
  }
}
/*--------------------------------------------------------------------------*/
int MeshUnitWithData::solveLinear(perulanganEnums::iterationstatus& status, const alat::GhostLinearSolver& B, const alat::GhostMatrix& A, alat::GhostVector& du, const alat::GhostVector& r) const
{
  getLinearSolver(B)->solve(status, A, du, r);
  return 1;
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::computeLinearSolver(perulanganEnums::matrixstatus& status, alat::GhostLinearSolver gB, alat::GhostMatrix gA, alat::GhostVector gu) const
{
  getLinearSolver(gB)->compute();
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::computeErrors(solvers::ErrorsMap& errormaps, std::string varname, const alat::GhostVector gu) const
{
  // std::cerr << "*** MeshUnitWithData::computeErrors()\n";
  int ivar = _var2index[varname];
  if(not _application->hasExactSolution(ivar)) return;
  const alat::VectorAllVariables& u = *getVector(gu);
  const FunctionInterface& exactsolutions = _application->getExactSolution(ivar);

  // std::cerr << "errormaps="<<errormaps;
  int ncomp = _vars[ivar].getNcomp();
  for(int iK=0; iK<_meshinfo->ncells;iK++)
  {
    _fems[ivar]->setVectorIndices(iK, _pdepartdata.vec_i[ivar], ncomp);
    u[ivar]->extract(_pdepartdata.vec_i[ivar], _pdepartdata.uloc[ivar]);
    _fems[ivar]->computeErrors(iK, errormaps, _pdepartdata.uloc[ivar], exactsolutions);
  }
  // std::cerr << "errormaps="<<errormaps;
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::preparePdePartsData(const alat::armaivec& ivars) const
{
  for(int ii=0; ii<ivars.size(); ii++)
  {
    int ivar = ivars[ii];
    assert(_fems[ivar]->getNPerCellConstant());
  }
  _nlocals.set_size(ivars.size());
  _ncomps.set_size(ivars.size());
  for(int ii=0; ii<ivars.size(); ii++)
  {
    int ivar = ivars[ii];
    _nlocals[ii] = _fems[ivar]->getNPerCell();
    _ncomps[ii] = _vars[ivar].getNcomp();
  }
  _pdepartdata.set_sizes(_ncomps, _nlocals);
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::computeRhs(alat::GhostVector gf, const alat::GhostVector gu) const
{
  // std::cerr << "*** MeshUnitWithData::computeRhs()\n";
  alat::VectorAllVariables& f = *getVector(gf);
  const alat::VectorAllVariables& u = *getVector(gu);
  int nvars = _vars.size();
  f.fillzeros();

  for(PdePartsMap::const_iterator p=_pdeparts.begin(); p!=_pdeparts.end();p++)
  {
    const alat::armaivec& ivars = p->second->getIvars();
    preparePdePartsData(ivars);
    if(p->second->loopCells())
    {
      for(int iK=0; iK<_meshinfo->ncells;iK++)
      {
        for(int ii=0; ii<ivars.size(); ii++)
        {
          int ivar = ivars[ii];
          int ncomp = _vars[ivar].getNcomp();
          _fems[ivar]->setVectorIndices(iK, _pdepartdata.vec_i[ii], ncomp);
          u[ivar]->extract(_pdepartdata.vec_i[ii], _pdepartdata.uloc[ii]);
          _fems[ivar]->setCell(iK);
          _pdepartdata.floc[ii].fill(arma::fill::zeros);
        }
        // for(int ii=0; ii<ivars.size(); ii++)
        // {
        //   int ivar = ivars[ii];
        //   std::cerr << "MeshUnitWithData::computeRhs() _pdepartdata.floc[ivar]" << _pdepartdata.floc[ivar];
        //   std::cerr << "MeshUnitWithData::computeRhs() _pdepartdata.uloc[ivar]" << _pdepartdata.uloc[ivar];
        // }
        p->second->computeRhsCell(iK, _pdepartdata.floc, _pdepartdata.uloc);
        for(int ii=0; ii<ivars.size(); ii++)
        {
          int ivar = ivars[ii];
          // std::cerr << "MeshUnitWithData::computeRhs() _pdepartdata.floc[ivar]" << _pdepartdata.floc[ivar];
          f[ivar]->assemble(_pdepartdata.vec_i[ii], _pdepartdata.floc[ii]);
        }
      }
    }
    // for(int ivar=0; ivar<nvars; ivar++)
    // {
    //   const alat::VectorOneVariable* fv = dynamic_cast<const alat::VectorOneVariable*>(f.get(ivar));
    //   std::cerr << "111 MeshUnitWithData::computeRhs() fv="<<*fv<<"\n";
    // }
    if(p->second->loopBoundary())
    {
      preparePdePartsData(ivars);
      for(mesh::MeshUnitInterface::BoundaryInformationMap::const_iterator pbdry=_meshinfo->bdrymesheunitsmap.begin();pbdry!=_meshinfo->bdrymesheunitsmap.end();pbdry++)
      {
        int color = pbdry->first;
        const alat::armaimat& cells_on_bdry = pbdry->second.getCellsOnBdryOfPlain();
        for(int i = 0; i < cells_on_bdry.n_cols; i++)
        {
          int iK = cells_on_bdry(0,i);
          int iS = cells_on_bdry(1,i);
          int iil = cells_on_bdry(2,i);
          for(int ii=0; ii<ivars.size(); ii++)
          {
            int ivar = ivars[ii];
            int ncomp = _vars[ivar].getNcomp();
            _fems[ivar]->setVectorIndices(iK, _pdepartdata.vec_i[ii], ncomp);
            u[ivar]->extract(_pdepartdata.vec_i[ii], _pdepartdata.uloc[ii]);
            _fems[ivar]->setCellBdry(iK, iS, iil);
            _pdepartdata.floc[ii].fill(arma::fill::zeros);
          }
          p->second->computeRhsBdry(color, iK, iS, iil, _pdepartdata.floc, _pdepartdata.uloc);
          for(int ii=0; ii<ivars.size(); ii++)
          {
            int ivar = ivars[ii];
            f[ivar]->assemble(_pdepartdata.vec_i[ii], _pdepartdata.floc[ii]);
          }
        }
      }
      // for(int ivar=0; ivar<nvars; ivar++)
      // {
      //   const alat::VectorOneVariable* fv = dynamic_cast<const alat::VectorOneVariable*>(f.get(ivar));
      //   std::cerr << "222 MeshUnitWithData::computeRhs() fv="<<*fv<<"\n";
      // }
    }
  }
  // {const alat::VectorOneVariable* fv = dynamic_cast<const alat::VectorOneVariable*>(f.get(0));
  // std::cerr << "MeshUnitWithData::computeRhs f=" << *fv << "\n";}
  // {const alat::VectorOneVariable* fv = dynamic_cast<const alat::VectorOneVariable*>(f.get(0));
  // std::cerr << "MeshUnitWithData::computeRhs f=" << *fv << "\n";}
  for(int ivar=0; ivar<nvars; ivar++)
  {
    if(_application->hasDirichlet(ivar))
    _fems[ivar]->strongDirichlet(f[ivar].get(), _application->getDirichlet(ivar), _application->getStrongDirichletColor());
    const alat::VectorOneVariable* fv = dynamic_cast<const alat::VectorOneVariable*>(f.get(ivar));
    if(not fv->is_finite())
    {
      std::cerr << "MeshUnitWithData::computeRhs() fv="<<*fv<<"\n";
      assert(0);
    }
  }
  // const alat::VectorOneVariable* fv = dynamic_cast<const alat::VectorOneVariable*>(f.get(0));
  // std::cerr << "MeshUnitWithData::computeRhs f=" << *fv << "\n";
}

/*--------------------------------------------------------------------------*/
void MeshUnitWithData::residual(perulanganEnums::residualstatus& status, alat::GhostVector& gr, const alat::GhostVector& gu, const alat::GhostVector& gf) const
{
  // std::cerr << "*** MeshUnitWithData::residual()\n";
  alat::VectorAllVariables& r = *getVector(gr);
  const alat::VectorAllVariables& u = *getVector(gu);
  const alat::VectorAllVariables& f = *getVector(gf);
  int nvars = _vars.size();
  assert(nvars==u.size());
  assert(nvars==f.size());
  assert(nvars==r.size());

  r.equal(&f);

  const alat::VectorOneVariable* rv = dynamic_cast<const alat::VectorOneVariable*>(r.get(0));
  const alat::VectorOneVariable* uv = dynamic_cast<const alat::VectorOneVariable*>(u.get(0));
  // std::cerr << "MeshUnitWithData::residual() f="<< *rv << "\n";
  // std::cerr << "MeshUnitWithData::residual() u="<< *uv << "\n";

  for(PdePartsMap::const_iterator p=_pdeparts.begin(); p!=_pdeparts.end();p++)
  {
    const alat::armaivec& ivars = p->second->getIvars();
    preparePdePartsData(ivars);
    // std::cerr << "MeshUnitWithData::residual() ivars="<<ivars.t();
    // Cells
    if(p->second->loopCells())
    {
      for(int iK=0; iK<_meshinfo->ncells;iK++)
      {
        // std::cerr << "*** MeshUnitWithData::residual() iK= << "<< iK << "\n";
        for(int ii=0; ii<ivars.size(); ii++)
        {
          int ivar = ivars[ii];
          int ncomp = _vars[ivar].getNcomp();
          _fems[ivar]->setVectorIndices(iK, _pdepartdata.vec_i[ii], ncomp);
          u[ivar]->extract(_pdepartdata.vec_i[ii], _pdepartdata.uloc[ii]);
          _fems[ivar]->setCell(iK);
          _pdepartdata.floc[ii].fill(arma::fill::zeros);
        }
        p->second->computeResidualCell(iK, _pdepartdata.floc, _pdepartdata.uloc);
        for(int ii=0; ii<ivars.size(); ii++)
        {
          int ivar = ivars[ii];
          std::string varname = _varnames[ivar];
          // std::cerr << "iK="<<iK<<" vec_i="<<_pdepartdata.vec_i[ivar].t()<<"  floc="<<_pdepartdata.floc[ivar];
          r[ivar]->assemble(_pdepartdata.vec_i[ii], _pdepartdata.floc[ii], -1.0);
        }
      }
    }
    // InteriorSides
    if(p->second->loopInteriorSides())
    {
      preparePdePartsData(ivars);
      // assert(0);
      for(int iS=0; iS<_meshinfo->nsides;iS++)
      {
        int iKin = _meshinfo->cells_of_sides(0,iS);
        int iKex = _meshinfo->cells_of_sides(1,iS);
        if(iKex<0) continue;
        // std::cerr << "*** MeshUnitWithData::residual() iKin, iKex= << "<< iKin << ","<<iKex<<"\n";
        if(not p->second->interiorsidecoupling(iKin, iKex)) continue;
        // std::cerr << "iS="<<iS<< " iKin="<<iKin<< " iKex="<<iKex << " nocoupling="<<nocoupling<<"\n";
        for(int ii=0; ii<ivars.size(); ii++)
        {
          int ivar = ivars[ii];
          int ncomp = _vars[ivar].getNcomp();
          _fems[ivar]->setVectorIndices(iKin, _pdepartdata.vec_i[ii], ncomp);
          assert(_femsex[ivar]);
          _femsex[ivar]->setVectorIndices(iKex, _pdepartdata.vec_iex[ii], ncomp);
          u[ivar]->extract(_pdepartdata.vec_i[ii], _pdepartdata.uloc[ii]);
          u[ivar]->extract(_pdepartdata.vec_iex[ii], _pdepartdata.ulocex[ii]);
          _fems[ivar]->setCell(iKin);
          _femsex[ivar]->setCell(iKex);
          _pdepartdata.floc[ii].fill(arma::fill::zeros);
          _pdepartdata.flocex[ii].fill(arma::fill::zeros);
        }
        p->second->computeResidualInteriorSide(iS, iKin, iKex, _pdepartdata.floc, _pdepartdata.flocex, _pdepartdata.uloc, _pdepartdata.ulocex);
        for(int ii=0; ii<ivars.size(); ii++)
        {
          int ivar = ivars[ii];
          std::string varname = _varnames[ivar];
          // std::cerr << "iK="<<iK<<" vec_i="<<_pdepartdata.vec_i[ivar].t()<<"  floc="<<_pdepartdata.floc[ivar];
          r[ivar]->assemble(_pdepartdata.vec_i[ii], _pdepartdata.floc[ii], -1.0);
          r[ivar]->assemble(_pdepartdata.vec_iex[ii], _pdepartdata.flocex[ii], -1.0);
        }
      }
    }
    // Boundary
    if(p->second->loopBoundary())
    {
      preparePdePartsData(ivars);
      for(mesh::MeshUnitInterface::BoundaryInformationMap::const_iterator pbdry=_meshinfo->bdrymesheunitsmap.begin();pbdry!=_meshinfo->bdrymesheunitsmap.end();pbdry++)
      {
        int color = pbdry->first;
        const alat::armaimat& cells_on_bdry = pbdry->second.getCellsOnBdryOfPlain();
        for(int i = 0; i < cells_on_bdry.n_cols; i++)
        {
          int iK = cells_on_bdry(0,i);
          int iS = cells_on_bdry(1,i);
          int iil = cells_on_bdry(2,i);
          for(int ii=0; ii<ivars.size(); ii++)
          {
            int ivar = ivars[ii];
            int ncomp = _vars[ivar].getNcomp();
            _fems[ivar]->setVectorIndices(iK, _pdepartdata.vec_i[ii], ncomp);
            u[ivar]->extract(_pdepartdata.vec_i[ivar], _pdepartdata.uloc[ii]);
            _fems[ivar]->setCellBdry(iK, iS, iil);
            _pdepartdata.floc[ii].fill(arma::fill::zeros);
          }
          p->second->computeResidualBdry(color, iK, iS, iil, _pdepartdata.floc, _pdepartdata.uloc);
          for(int ii=0; ii<ivars.size(); ii++)
          {
            int ivar = ivars[ii];
            r[ivar]->assemble(_pdepartdata.vec_i[ii], _pdepartdata.floc[ii], -1.0);
          }
        }
      }
    }
    // Global
    p->second->computeResidualGlobal(r, u);
  }
  for(int ivar=0; ivar<nvars; ivar++)
  {
    _fems[ivar]->strongDirichletZero(r[ivar].get(), _application->getStrongDirichletColor());
  }
  // std::cerr << "MeshUnitWithData::residual() r="<< *rv << "\n";
}

// /*--------------------------------------------------------------------------*/
// void MeshUnitWithData::prepareCellMatrix(int iK, const alat::armaivec& ivars) const
// {
//   _pdepartdata.aloc.set_size(ivars.size(),ivars.size());
//   // _pdepartdata.aloc_i.set_size(ivars.size(),ivars.size());
//   // _pdepartdata.aloc_j.set_size(ivars.size(),ivars.size());
//   for(int ii=0; ii<ivars.size(); ii++)
//   {
//     int ivar = ivars[ii];
//     int ncompi = _fems[ivar]->getNcomp();
//     int nglobi = _fems[ivar]->getN();
//     int nloccelli = _fems[ivar]->getNPerCell(iK);
//     alat::armaivec indicesi(nloccelli);
//     _fems[ivar]->indicesOfCell(iK, indicesi);
//     for(int jj=0; jj<ivars.size(); jj++)
//     {
//       int jvar = ivars[jj];
//       int ncompj = _fems[jvar]->getNcomp();
//       int nglobj = _fems[jvar]->getN();
//       int nloccellj = _fems[jvar]->getNPerCell(iK);
//       alat::armaivec indicesj(nloccellj);
//       _fems[jvar]->indicesOfCell(iK, indicesj);
//       _pdepartdata.aloc(ii,jj).set_size(ncompi*ncompj*nloccelli*nloccellj);
//       _pdepartdata.aloc(ii,jj).fill(arma::fill::zeros);
//       int count=0;
//       for(int iloc=0; iloc<nloccelli;iloc++)
//       {
//         int i = indicesi[iloc];
//         for(int jloc=0; jloc<nloccellj;jloc++)
//         {
//           int j = indicesj[jloc];
//           for(int icomp=0;icomp<ncompi;icomp++)
//           {
//             for(int jcomp=0;jcomp<ncompj;jcomp++)
//             {
//               _pdepartdata.aloc_i(ii,jj)[count] = icomp*nglobi + i;
//               _pdepartdata.aloc_j(ii,jj)[count] = jcomp*nglobj + j;
//               count++;
//             }
//           }
//         }
//       }
//     }
//   }
// }
//
// /*--------------------------------------------------------------------------*/
// void MeshUnitWithData::prepareCellMatrix(int iKin, int iKex, const alat::armaivec& ivars) const
// {
//   _pdepartdata.aloc.set_size(ivars.size(),ivars.size());
//   _pdepartdata.aloc_inex.set_size(ivars.size(),ivars.size());
//   _pdepartdata.aloc_exin.set_size(ivars.size(),ivars.size());
//   _pdepartdata.aloc_exex.set_size(ivars.size(),ivars.size());
//   _pdepartdata.aloc_i.set_size(ivars.size(),ivars.size());
//   _pdepartdata.aloc_j.set_size(ivars.size(),ivars.size());
//   _pdepartdata.aloc_i_ex.set_size(ivars.size(),ivars.size());
//   _pdepartdata.aloc_j_ex.set_size(ivars.size(),ivars.size());
//   for(int ii=0; ii<ivars.size(); ii++)
//   {
//     int ivar = ivars[ii];
//     int ncompi = _fems[ivar]->getNcomp();
//     int nglobi = _fems[ivar]->getN();
//     int nloccelli = _fems[ivar]->getNPerCell(iKin);
//     alat::armaivec indicesi_in(nloccelli), indicesi_ex(nloccelli);
//     _fems[ivar]->indicesOfCell(iKin, indicesi_in);
//     _femsex[ivar]->indicesOfCell(iKex, indicesi_ex);
//     for(int jj=0; jj<ivars.size(); jj++)
//     {
//       int jvar = ivars[jj];
//       int ncompj = _fems[jvar]->getNcomp();
//       int nglobj = _fems[jvar]->getN();
//       int nloccellj = _fems[jvar]->getNPerCell(iKin);
//       alat::armaivec indicesj_in(nloccellj), indicesj_ex(nloccellj);
//       _fems[jvar]->indicesOfCell(iKin, indicesj_in);
//       _femsex[jvar]->indicesOfCell(iKex, indicesj_ex);
//       _pdepartdata.aloc     (ii,jj).set_size(ncompi*ncompj*nloccelli*nloccellj);
//       _pdepartdata.aloc     (ii,jj).fill(arma::fill::zeros);
//       _pdepartdata.aloc_inex(ii,jj).set_size(ncompi*ncompj*nloccelli*nloccellj);
//       _pdepartdata.aloc_inex(ii,jj).fill(arma::fill::zeros);
//       _pdepartdata.aloc_exin(ii,jj).set_size(ncompi*ncompj*nloccelli*nloccellj);
//       _pdepartdata.aloc_exin(ii,jj).fill(arma::fill::zeros);
//       _pdepartdata.aloc_exex(ii,jj).set_size(ncompi*ncompj*nloccelli*nloccellj);
//       _pdepartdata.aloc_exex(ii,jj).fill(arma::fill::zeros);
//       int count=0;
//       for(int iloc=0; iloc<nloccelli;iloc++)
//       {
//         int iin = indicesi_in[iloc];
//         int iex = indicesi_ex[iloc];
//         for(int jloc=0; jloc<nloccellj;jloc++)
//         {
//           int jin = indicesj_in[jloc];
//           int jex = indicesj_ex[jloc];
//           for(int icomp=0;icomp<ncompi;icomp++)
//           {
//             for(int jcomp=0;jcomp<ncompj;jcomp++)
//             {
//               _pdepartdata.aloc_i   (ii, jj)[count] = icomp*nglobi + iin;
//               _pdepartdata.aloc_j   (ii, jj)[count] = jcomp*nglobj + jin;
//               _pdepartdata.aloc_i_ex(ii, jj)[count] = icomp*nglobi + iex;
//               _pdepartdata.aloc_j_ex(ii, jj)[count] = jcomp*nglobj + jex;
//               count++;
//             }
//           }
//         }
//       }
//     }
//   }
// }
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::computeJacobian(perulanganEnums::matrixstatus& status, alat::GhostMatrix gA, const alat::GhostVector gu) const
{
  // std::cerr << "*** MeshUnitWithData::computeJacobian()\n";
  alat::MatrixAllVariables& A = *getMatrix(gA);
  const alat::VectorAllVariables& u = *getVector(gu);
  int nvars = _vars.size();
  A.fillzeros();

  // const alat::MatrixOneVariable* Av = dynamic_cast<const alat::MatrixOneVariable*>(A(0,0).get());
  // std::cerr << "MeshUnitWithData::computeJacobian A=" << *Av << "\n";

  for(PdePartsMap::const_iterator p=_pdeparts.begin(); p!=_pdeparts.end();p++)
  {
    const alat::armaivec& ivars = p->second->getIvars();
    preparePdePartsData(ivars);
    // Cells
    if(p->second->loopCells())
    {
      for(int iK=0; iK<_meshinfo->ncells;iK++)
      {
        for(int ii=0; ii<ivars.size(); ii++)
        {
          int ivar = ivars[ii];
          int ncomp = _vars[ivar].getNcomp();
          _fems[ivar]->setVectorIndices(iK, _pdepartdata.vec_i[ii], ncomp);
          u[ivar]->extract(_pdepartdata.vec_i[ii], _pdepartdata.uloc[ii]);
          _fems[ivar]->setCell(iK);
          for(int jj=0; jj<ivars.size(); jj++)
          {
            _pdepartdata.aloc(ii,jj).fill(arma::fill::zeros);
          }
        }
        p->second->computeMatrixCell(iK, _pdepartdata.aloc,_pdepartdata.uloc);
        for(int ii=0; ii<ivars.size(); ii++)
        {
          int ivar = ivars[ii];
          for(int jj=0; jj<ivars.size(); jj++)
          {
            int jvar = ivars[jj];
            // std::cerr << "_pdepartdata.vec_i[ii]="<<_pdepartdata.vec_i[ii];
            // std::cerr << "_pdepartdata.aloc(ii,jj)="<<_pdepartdata.aloc(ii,jj);
            A(ivar,jvar)->assemble(_pdepartdata.aloc(ii,jj), _pdepartdata.vec_i[ii], _pdepartdata.vec_i[jj]);
          }
        }
      }
    }
    // InteriorSides
    if(p->second->loopInteriorSides())
    {
      preparePdePartsData(ivars);
      for(int iS=0; iS<_meshinfo->nsides;iS++)
      {
        int iKin = _meshinfo->cells_of_sides(0,iS);
        int iKex = _meshinfo->cells_of_sides(1,iS);
        if(iKex<0) continue;
        if(not p->second->interiorsidecoupling(iKin, iKex)) continue;
        for(int ii=0; ii<ivars.size(); ii++)
        {
          int ivar = ivars[ii];
          int ncomp = _vars[ivar].getNcomp();
          _fems[ivar]->setVectorIndices(iKin, _pdepartdata.vec_i[ii], ncomp);
          _femsex[ivar]->setVectorIndices(iKex, _pdepartdata.vec_iex[ii], ncomp);
          u[ivar]->extract(_pdepartdata.vec_i[ii], _pdepartdata.uloc[ii]);
          u[ivar]->extract(_pdepartdata.vec_iex[ii], _pdepartdata.ulocex[ii]);
          _fems[ivar]->setCell(iKin);
          _femsex[ivar]->setCell(iKex);
          _pdepartdata.floc[ii].fill(arma::fill::zeros);
          _pdepartdata.flocex[ii].fill(arma::fill::zeros);
          for(int jj=0; jj<ivars.size(); jj++)
          {
            _pdepartdata.aloc     (ii,jj).fill(arma::fill::zeros);
            _pdepartdata.aloc_inex(ii,jj).fill(arma::fill::zeros);
            _pdepartdata.aloc_exin(ii,jj).fill(arma::fill::zeros);
            _pdepartdata.aloc_exex(ii,jj).fill(arma::fill::zeros);
          }
        }
        p->second->computeMatrixInteriorSide(iS, iKin, iKex, _pdepartdata.aloc, _pdepartdata.aloc_inex, _pdepartdata.aloc_exin, _pdepartdata.aloc_exex,_pdepartdata.uloc, _pdepartdata.ulocex);
        for(int ii=0; ii<ivars.size(); ii++)
        {
          int ivar = ivars[ii];
          for(int jj=0; jj<ivars.size(); jj++)
          {
            int jvar = ivars[jj];
            A(ivar,jvar)->assemble(_pdepartdata.aloc     (ii,jj), _pdepartdata.vec_i[ii]  , _pdepartdata.vec_i[jj]);
            A(ivar,jvar)->assemble(_pdepartdata.aloc_inex(ii,jj), _pdepartdata.vec_i[ii]  , _pdepartdata.vec_iex[jj]);
            A(ivar,jvar)->assemble(_pdepartdata.aloc_exin(ii,jj), _pdepartdata.vec_iex[ii], _pdepartdata.vec_i[jj]);
            A(ivar,jvar)->assemble(_pdepartdata.aloc_exex(ii,jj), _pdepartdata.vec_iex[ii], _pdepartdata.vec_iex[jj]);
          }
        }
      }
    }
    // Boundary
    if(p->second->loopBoundary())
    {
      preparePdePartsData(ivars);
      for(mesh::MeshUnitInterface::BoundaryInformationMap::const_iterator pbdry=_meshinfo->bdrymesheunitsmap.begin();pbdry!=_meshinfo->bdrymesheunitsmap.end();pbdry++)
      {
        int color = pbdry->first;
        const alat::armaimat& cells_on_bdry = pbdry->second.getCellsOnBdryOfPlain();
        for(int i = 0; i < cells_on_bdry.n_cols; i++)
        {
          int iK = cells_on_bdry(0,i);
          int iS = cells_on_bdry(1,i);
          int iil = cells_on_bdry(2,i);
          for(int ii=0; ii<ivars.size(); ii++)
          {
            int ivar = ivars[ii];
            int ncomp = _vars[ivar].getNcomp();
            _fems[ivar]->setVectorIndices(iK, _pdepartdata.vec_i[ivar], ncomp);
            u[ivar]->extract(_pdepartdata.vec_i[ivar], _pdepartdata.uloc[ivar]);
            _fems[ivar]->setCellBdry(iK, iS, iil);
            for(int jj=0; jj<ivars.size(); jj++)
            {
              _pdepartdata.aloc(ii,jj).fill(arma::fill::zeros);
            }
          }
          p->second->computeMatrixBdry(color, iK, iS, iil, _pdepartdata.aloc,_pdepartdata.uloc);
          for(int ii=0; ii<ivars.size(); ii++)
          {
            int ivar = ivars[ii];
            for(int jj=0; jj<ivars.size(); jj++)
            {
              int jvar = ivars[jj];
              A(ivar,jvar)->assemble(_pdepartdata.aloc(ii,jj), _pdepartdata.vec_i[ii], _pdepartdata.vec_i[jj]);
            }
          }
        }
      }
    }
    // Global
    p->second->computeMatrixGlobal(A, u);
  }
  for(int ivar=0; ivar<nvars; ivar++)
  {
    for(int jvar=0; jvar<nvars; jvar++)
    {
      if(ivar==jvar)
      {
        _fems[ivar]->strongDirichlet(ivar, A, _application->getStrongDirichletColor());
      }
    }
  }
  // for(int ivar=0; ivar<nvars; ivar++)
  // {
  //   for(int jvar=0; jvar<nvars; jvar++)
  //   {
  //     const alat::MatrixOneVariable* Av = dynamic_cast<const alat::MatrixOneVariable*>(A(ivar,jvar).get());
  //     std::cerr << "MeshUnitWithData::computeJacobian A= ivar/jvar" <<ivar << "/" << jvar << " "<< *Av << "\n";
  //   }
  // }
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::vectorZero(alat::GhostVector& gu) const
{
  getVector(gu)->fillzeros();
}
void MeshUnitWithData::vectorAdd(alat::GhostVector& gu, double s, const alat::GhostVector& gv) const
{
  getVector(gu)->add(s, getVector(gv));
}
void MeshUnitWithData::vectorEqual(alat::GhostVector& gu, const alat::GhostVector& gv) const
{
  getVector(gu)->equal(getVector(gv));
}
double MeshUnitWithData::vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const
{
  return getVector(gu)->dot(getVector(gv));
}

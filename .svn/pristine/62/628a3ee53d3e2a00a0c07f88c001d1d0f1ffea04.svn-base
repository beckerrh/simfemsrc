#include  "Alat/stringvector.hpp"
#include  "Solvers/applicationinterface.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
ApplicationInterface::~ApplicationInterface() {}
ApplicationInterface::ApplicationInterface(): alat::InterfaceBase(), _mesh(NULL){}
ApplicationInterface::ApplicationInterface( const ApplicationInterface& applicationinterface): alat::InterfaceBase(applicationinterface)
{
  assert(0);
}
ApplicationInterface& ApplicationInterface::operator=( const ApplicationInterface& applicationinterface)
{
  assert(0);
  alat::InterfaceBase::operator=(applicationinterface);
  return *this;
}
std::string ApplicationInterface::getClassName() const
{
  return "ApplicationInterface";
}
void ApplicationInterface::setModel(const solvers::ModelInterface* model){_model=model;}

/*--------------------------------------------------------------------------*/
std::unique_ptr<solvers::InitialConditionInterface> ApplicationInterface::newInitialCondition(std::string varname) const
{
  return std::unique_ptr<solvers::InitialConditionInterface>(nullptr);
}
std::unique_ptr<solvers::RightHandSideInterface> ApplicationInterface::newRightHandSide(std::string varname) const
{
  return std::unique_ptr<solvers::RightHandSideInterface>(nullptr);
}
std::unique_ptr<solvers::DirichletInterface> ApplicationInterface::newDirichlet(std::string varname) const
{
  return std::unique_ptr<solvers::DirichletInterface>(nullptr);
}
std::unique_ptr<solvers::NeumannInterface> ApplicationInterface::newNeumann(std::string varname) const
{
  return std::unique_ptr<solvers::NeumannInterface>(nullptr);
}
std::shared_ptr<solvers::FunctionInterface> ApplicationInterface::newDataFunction(std::string varname) const
{
  return std::shared_ptr<solvers::FunctionInterface>(nullptr);
}
std::shared_ptr<solvers::FunctionInterface> ApplicationInterface::newExactSolution(std::string varname) const
{
  return std::shared_ptr<solvers::FunctionInterface>(nullptr);
}

const InitialConditionInterface& ApplicationInterface::getInitialCondition(int ivar)const{return *_initialconditions[ivar].get();}
const RightHandSideInterface& ApplicationInterface::getRightHandSide(int ivar)const{assert(ivar<_righthandsides.size());return *_righthandsides[ivar].get();}
const DirichletInterface& ApplicationInterface::getDirichlet(int ivar)const{assert(ivar<_dirichlets.size()); return *_dirichlets[ivar].get();}
const NeumannInterface& ApplicationInterface::getNeumann(int ivar)const{assert(ivar<_neumanns.size()); return *_neumanns[ivar].get();}
const FunctionInterface& ApplicationInterface::getDataFunction(int ivar)const{assert(ivar<_datafunctions.size()); return *_datafunctions[ivar].get();}
const solvers::FunctionInterface& ApplicationInterface::getExactSolution(int ivar)const
{
  assert(_exactsolutions[ivar]);
  return *_exactsolutions[ivar].get();
}
solvers::FunctionInterface& ApplicationInterface::getExactSolution(int ivar)
{
  assert(_exactsolutions[ivar]);
  return *_exactsolutions[ivar].get();
}

const InitialConditionInterface& ApplicationInterface::getInitialCondition(std::string varname)const
{
  return *_initialconditions[_indexofvar[varname]].get();
}
const RightHandSideInterface& ApplicationInterface::getRightHandSide(std::string varname)const
{
  return *_righthandsides[_indexofvar[varname]].get();
}
const DirichletInterface& ApplicationInterface::getDirichlet(std::string varname)const
{
  return *_dirichlets[_indexofvar[varname]].get();
}
const NeumannInterface& ApplicationInterface::getNeumann(std::string varname)const
{
  return *_neumanns[_indexofvar[varname]].get();
}
const FunctionInterface& ApplicationInterface::getDataFunction(std::string varname)const
{
  assert(_datafunctions[_indexofdatavar[varname]]);
  return *_datafunctions[_indexofdatavar[varname]];
}
const solvers::FunctionInterface& ApplicationInterface::getExactSolution(std::string varname)const
{
  return *_exactsolutions[_indexofvar[varname]].get();
}
solvers::FunctionInterface& ApplicationInterface::getExactSolution(std::string varname)
{
  return *_exactsolutions[_indexofvar[varname]].get();
}


bool ApplicationInterface::hasInitialCondition(int ivar)const
{
  assert(ivar<_initialconditions.size());
  return _initialconditions[ivar] != nullptr;
}
bool ApplicationInterface::hasDataFunction(int ivar)const
{
  assert(ivar<_datafunctions.size());
  return _datafunctions[ivar] != nullptr;
}
bool ApplicationInterface::hasExactSolution(int ivar)const
{
  assert(ivar<_exactsolutions.size());
  return _exactsolutions[ivar] != nullptr;
}
bool ApplicationInterface::hasExactSolution(std::string varname)const
{
  return hasExactSolution(_indexofvar[varname]);
}

bool ApplicationInterface::hasRightHandSide(int ivar)const{return _righthandsides[ivar] != nullptr;}
bool ApplicationInterface::hasDirichlet(int ivar)const{return _dirichlets[ivar] != nullptr;}
const alat::IntSet& ApplicationInterface::getStrongDirichletColor() const {return _dircolors;}

/*--------------------------------------------------------------------------*/
void ApplicationInterface::initApplication(const mesh::MeshUnitInterface* mesh, const alat::StringVector& varnames, const alat::StringVector& varnamesdata, const alat::armaivec& ncomps, const alat::armaivec& ncompsdata)
{
  _mesh=mesh;
  int nvars = varnames.size();
  _initialconditions.set_size(nvars);
  _righthandsides.set_size(nvars);
  _dirichlets.set_size(nvars);
  _neumanns.set_size(nvars);
  _exactsolutions.set_size(nvars);
  int dim = _mesh->getDimension();
  for(int ivar=0;ivar<nvars;ivar++)
  {
    std::string varname = varnames[ivar];
    _indexofvar[varname] = ivar;
    _exactsolutions[ivar] = newExactSolution(varname);
    _initialconditions[ivar] = newInitialCondition(varname);
    _righthandsides[ivar] = newRightHandSide(varname);
    _dirichlets[ivar] = newDirichlet(varname);
    _neumanns[ivar] = newNeumann(varname);
    if(_initialconditions[ivar]) {_initialconditions[ivar]->setNcompAndDim(ncomps[ivar], dim);}
    if(_righthandsides[ivar]) {_righthandsides[ivar]->setNcompAndDim(ncomps[ivar], dim);}
    if(_dirichlets[ivar]) {_dirichlets[ivar]->setNcompAndDim(ncomps[ivar], dim);}
    if(_neumanns[ivar]) {_neumanns[ivar]->setNcompAndDim(ncomps[ivar], dim);}
  }
  int nvarsdata = varnamesdata.size();
  _datafunctions.set_size(nvarsdata);
  for(int ivar=0;ivar<nvarsdata;ivar++)
  {
    std::string varname = varnamesdata[ivar];
    _indexofdatavar[varname] = ivar;
    _datafunctions[ivar] = newDataFunction(varname);
    if(_datafunctions[ivar]) {_datafunctions[ivar]->setNcompAndDim(ncompsdata[ivar], dim);}
  }
  _dircolors.clear();
  for(mesh::MeshUnitInterface::BoundaryInformationMap::const_iterator p=_mesh->getBoundaryInformationMap().begin();p!=_mesh->getBoundaryInformationMap().end();p++)
  {
    int color = p->first;
    if(isStrongDirichlet(color)) _dircolors.insert(color);
  }
}
/*--------------------------------------------------------------------------*/
std::string ApplicationInterface::getInfo() const
{
  std::stringstream ss;
  ss << "\t RightHandSides: ";
  int nvars = _righthandsides.size();
  for(int ivar=0;ivar<nvars;ivar++)
  {
    if(_righthandsides[ivar])
    {
      ss << _righthandsides[ivar]->getClassName() << " ";
    }
    else
    {
      ss << " NULL ";
    }
  }
  ss << "\n\t Dirichlets: ";
  for(int ivar=0;ivar<nvars;ivar++)
  {
    if(_dirichlets[ivar])
    {
      ss << _dirichlets[ivar]->getClassName() << " ";
    }
    else
    {
      ss << " NULL ";
    }
  }
  ss << "\n\t Neumanns: ";
  for(int ivar=0;ivar<nvars;ivar++)
  {
    if(_neumanns[ivar])
    {
      ss << _neumanns[ivar]->getClassName() << " ";
    }
    else
    {
      ss << " NULL ";
    }
  }
  ss << "\n\t Exactsolutoins: ";
  for(int ivar=0;ivar<nvars;ivar++)
  {
    if(_exactsolutions[ivar])
    {
      ss << _exactsolutions[ivar]->getClassName() << " ";
    }
    else
    {
      ss << " NULL ";
    }
  }
  return ss.str();
}

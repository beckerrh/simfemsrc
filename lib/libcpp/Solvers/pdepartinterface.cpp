#include  "Alat/matrixonevariableinterface.hpp"
#include  "Solvers/meshunitwithdatainterface.hpp"
#include  "Solvers/pdepartinterface.hpp"
#include  "Solvers/variable.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
void PdePartData::set_size(const alat::armaivec& ncomps, const alat::armaivec& nlocals)
{
  int nvars = ncomps.size();
  uloc  .set_size(nvars);
  floc  .set_size(nvars);
  vec_i .set_size(nvars);
  aloc  .set_size(nvars, nvars);
  for(int ivar=0;ivar<nvars;ivar++)
  {
    int ncompi = ncomps[ivar];
    int nlocali = nlocals[ivar];
    int sizei = ncompi*nlocali;
    uloc [ivar].set_size(sizei);
    floc [ivar].set_size(sizei);
    vec_i[ivar].set_size(sizei);
    for(int jvar=0;jvar<nvars;jvar++)
    {
      int ncompj = ncomps[jvar];
      int nlocalj = nlocals[jvar];
      int sizej = ncompj*nlocalj;
      aloc  (ivar, jvar).set_size(sizei, sizej);
    }
  }
}

/*--------------------------------------------------------------------------*/
void PdePartData::set_sizes(const alat::armaivec& ncomps, const alat::armaivec& nlocals)
{
  int nvars = ncomps.size();
  uloc  .set_size(nvars);
  ulocex.set_size(nvars);
  floc  .set_size(nvars);
  flocex.set_size(nvars);
  vec_i .set_size(nvars);
  vec_iex.set_size(nvars);
  aloc     .set_size(nvars, nvars);
  aloc_inex.set_size(nvars, nvars);
  aloc_exin.set_size(nvars, nvars);
  aloc_exex.set_size(nvars, nvars);
  for(int ivar=0;ivar<nvars;ivar++)
  {
    int ncompi = ncomps[ivar];
    int nlocali = nlocals[ivar];
    int sizei = ncompi*nlocali;
    uloc  [ivar].set_size(sizei);
    ulocex[ivar].set_size(sizei);
    floc  [ivar].set_size(sizei);
    flocex[ivar].set_size(sizei);
    vec_i  [ivar].set_size(sizei);
    vec_iex[ivar].set_size(sizei);
    for(int jvar=0;jvar<nvars;jvar++)
    {
      int ncompj = ncomps[jvar];
      int nlocalj = nlocals[jvar];
      int sizej = ncompj*nlocalj;
      aloc     (ivar, jvar).set_size(sizei, sizej);
      aloc_inex(ivar, jvar).set_size(sizei, sizej);
      aloc_exin(ivar, jvar).set_size(sizei, sizej);
      aloc_exex(ivar, jvar).set_size(sizei, sizej);
    }
  }
}

/*--------------------------------------------------------------------------*/
PdePartInterface::~PdePartInterface() {}
PdePartInterface::PdePartInterface(alat::StringList vars, solver_options::pdepart::opts opts): alat::InterfaceBase(), _opts(opts), _vars(vars), _application(NULL), _model(NULL){}
PdePartInterface::PdePartInterface( const PdePartInterface& pdepartinterface): alat::InterfaceBase(pdepartinterface), _opts(pdepartinterface._opts)
{
  assert(0);
}
PdePartInterface& PdePartInterface::operator=( const PdePartInterface& pdepartinterface)
{
  assert(0);
  alat::InterfaceBase::operator=(pdepartinterface);
  return *this;
}
std::string PdePartInterface::getClassName() const
{
  return "PdePartInterface";
}
/*--------------------------------------------------------------------------*/
std::string PdePartInterface::getInfo() const
{
  std::stringstream ss;
  ss << getClassName() << " opts="<<_opts;
  return ss.str();
}
/*--------------------------------------------------------------------------*/
const alat::armaivec& PdePartInterface::getIvars() const {return _ivars;}

/*--------------------------------------------------------------------------*/
void PdePartInterface::setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters)
{
  // std::cerr << "_vars="<<_vars<<"\n";
  // std::cerr << "var2index="<<var2index<<"\n";
  _ivars.set_size(_vars.size());
  int count=0;
  for(alat::StringList::const_iterator p=_vars.begin();p!=_vars.end();p++)
  {
    _ivars[count++] = var2index[*p];
  }
  // std::cerr << "_ivars="<<_ivars<<"\n";
}

/*--------------------------------------------------------------------------*/
bool PdePartInterface::loopCells() const
{
  // std::cerr << "PdePartInterface::loopCells() _opts="<<_opts<<"\n";
  return solver_options::pdepart::has_cell(_opts);
}
bool PdePartInterface::loopBoundary() const
{
  return solver_options::pdepart::has_bdry(_opts);
}
bool PdePartInterface::loopInteriorSides() const
{
  return solver_options::pdepart::has_iside(_opts);
}
bool PdePartInterface::interiorsidecoupling(int iKin, int iKex) const
{
  return true;
}
/*--------------------------------------------------------------------------*/

void PdePartInterface::computeResidualInteriorSide(int iS, int iKin, int iKex, solvers::PdePartData::vec& flocin, solvers::PdePartData::vec& flocex, const solvers::PdePartData::vec& ulocin, const solvers::PdePartData::vec& ulocex){}
void PdePartInterface::computeMatrixInteriorSide(int iS, int iKin, int iKex, solvers::PdePartData::mat& matinin, solvers::PdePartData::mat& matinex, solvers::PdePartData::mat& matexin, solvers::PdePartData::mat& matexex, const solvers::PdePartData::vec& ulocin, const solvers::PdePartData::vec& ulocex)const{}

/*--------------------------------------------------------------------------*/
void PdePartInterface::rhsCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const
{_notWritten("rhsCell");}
void PdePartInterface::residualCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const
{_notWritten("residualCell");}
void PdePartInterface::matrixCell(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const
{_notWritten("matrixCell");}
void PdePartInterface::rhsBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const{_notWritten("rhsBdry");}
void PdePartInterface::residualBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const{_notWritten("residualBdry");}
void PdePartInterface::matrixBdry(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const{_notWritten("matrixBdry");}
void PdePartInterface::rhsBdryCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const
{
  rhsCell(floc, fems);
}
void PdePartInterface::residualBdryCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const
{
  residualCell(floc, fems);
}
void PdePartInterface::matrixBdryCell(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const
{
  matrixCell(mat, fems);
}
void PdePartInterface::prepareRhsCellBdry(int iK) const{}

void PdePartInterface::additionCouplings(alat::Matrix<alat::SparsityPatternSoft>& sparsitypatternsoft)const {}
void PdePartInterface::computeMatrixGlobal(alat::MatrixAllVariables& A, const alat::VectorAllVariables& u)const {}
void PdePartInterface::computeResidualGlobal(alat::VectorAllVariables& r, const alat::VectorAllVariables& u)const {}

/*--------------------------------------------------------------------------*/
void PdePartInterface::setTimeMeshData(double time, double dt) {_time=time; _dt=dt;}
/*--------------------------------------------------------------------------*/
void PdePartInterface::initPdePart(const solvers::MeshInfo* meshinfo, const alat::Map<std::string, int>& var2index, const solvers::MeshUnitWithDataInterface* meshunitwithdata, const solvers::Parameters& parameters)
{
  _mesh = meshunitwithdata->getMesh();
  _meshunit = meshunitwithdata;
  _meshinfo = meshinfo;
  _application = meshunitwithdata->getApplication();
  _model = meshunitwithdata->getModel();
  _fems  = meshunitwithdata->getFemMap() ;
  _femsex  = meshunitwithdata->getFemExMap() ;
  setData(var2index, parameters);
  _opts = setOptions();

  int maxpoints=0;
  int nvars = _fems->size();
  _femdatas.set_size(nvars);
  _femdatasex.set_size(nvars);
  for(int ivar=0; ivar<nvars; ivar++)
  {
    if((*_fems)[ivar]->noIntegration()) continue;
    int npoints = (*_fems)[ivar]->getFormula()->n();
    if(npoints>=maxpoints)
    {
      maxpoints=npoints;
      _femforintegration = (*_fems)[ivar].get();
    }
    _femdatas[ivar] = NULL;
    _femdatasex[ivar] = NULL;
  }

  _cellisbdry.set_size(meshinfo->ncells);
  _cellisbdry.fill(arma::fill::zeros);
  for(int ivar=0; ivar<nvars; ivar++)
  {
    (*_fems)[ivar]->setCellIsBdry(_cellisbdry);
  }
}

/*--------------------------------------------------------------------------*/
void PdePartInterface::computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const{_notWritten("computeRhsCell");}
void PdePartInterface::computeRhsBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const{_notWritten("computeRhsBdry");}
void PdePartInterface::computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const{_notWritten("computeResidualCell");}
void PdePartInterface::computeResidualBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const{_notWritten("computeResidualBdry");}
void PdePartInterface::computeMatrixCell(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const{_notWritten("computeMatrixCell");}
void PdePartInterface::computeMatrixBdry(int color, int iK, int iS, int iil, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const{_notWritten("computeMatrixBdry");}

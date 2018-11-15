#include  "Solvers/feminterface.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
FemInterface::~FemInterface() {}
FemInterface::FemInterface(): alat::InterfaceBase(), _mesh(NULL), _meshinfo(NULL){}
FemInterface::FemInterface( const FemInterface& feminterface): alat::InterfaceBase(feminterface)
{
  (*this).operator=(feminterface);
  _mesh = feminterface._mesh;
  _meshinfo = feminterface._meshinfo;
}
FemInterface& FemInterface::operator=( const FemInterface& feminterface)
{
  assert(0);
  alat::InterfaceBase::operator=(feminterface);
  return *this;
}
std::string FemInterface::getClassName() const
{
  return "FemInterface";
}
std::string FemInterface::getInfo() const
{
  std::stringstream ss;
  ss << getClassName();
  return ss.str();
}

/*--------------------------------------------------------------------------*/
void FemInterface::initFem(int ivar, const mesh::MeshUnitInterface* mesh, const MeshInfo* meshinfo, int ncomp)
{
  _ivar = ivar;
  _ncomp = ncomp;
  _mesh = mesh;
  _dim = mesh->getDimension();
  _meshinfo = meshinfo;
  initData();
}
int FemInterface::getNcomp() const {return _ncomp;}

/*--------------------------------------------------------------------------*/
int FemInterface::getNPerCell(int iK) const{_notWritten("getNPerCell");}
void FemInterface::indicesOfCell(int iK, alat::armaivec& indices) const{_notWritten("indicesOfCell");}
const solvers::IntegrationFormulaInterface* FemInterface::getFormula() const{_notWritten("getFormula");}
const solvers::IntegrationFormulaInterface* FemInterface::getFormulaErrors() const{_notWritten("getFormulaErrors");}
const solvers::IntegrationFormulaInterface* FemInterface::getFormulaBdry() const{_notWritten("getFormulaBdry");}
const solvers::IntegrationFormulaInterface* FemInterface::getFormulaRhs() const{_notWritten("getFormulaRhs");}
void FemInterface::setCell(int iK){_notWritten("setCell");}
void FemInterface::setCellBdry(int iK, int iS, int iil){_notWritten("setCellIsBdry");}
const FemData& FemInterface::referencePoint(const alat::Node& vhat, double weight){_notWritten("referencePoint");}
const FemData& FemInterface::referencePointWithData(const alat::Node& vhat, double weight, const arma::mat& uloc){_notWritten("referencePointWithData");}
const FemData& FemInterface::referencePointBdry(const alat::Node& vhat, double weight){_notWritten("referencePointBdry");}
const FemData& FemInterface::referencePointBdryWithData(const alat::Node& vhat, double weight, const arma::mat& uloc){_notWritten("referencePointBdryWithData");}
const FemData& FemInterface::referencePointBdryCellWithData(const alat::Node& vhat, double weight, const arma::mat& uloc){_notWritten("referencePointBdryCellWithData");}
void FemInterface::setVectorIndices(int iK, alat::armaimat& vec_i)const{_notWritten("setVectorIndices");}

void FemInterface::computeGrad(arma::mat& ugrad, const arma::mat& uloc) const{_notWritten("computeGrad");}
void FemInterface::computeFunction(arma::vec& u, const arma::mat& uloc) const{_notWritten("computeFunction");}
const solvers::FemData& FemInterface::getFemdata() const{_notWritten("getFemdata");}


bool FemInterface::canInterpolateToP1()const{return true;}
void FemInterface::setCellIsBdry(arma::uvec& cellisbdry){}
void FemInterface::setIsi(int iK){_notWritten("setIsi");}
void FemInterface::strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const{_notWritten("strongDirichlet");}
void FemInterface::strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const{_notWritten("strongDirichletZero");}
void FemInterface::strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const{_notWritten("strongDirichlet");}
void FemInterface::interpolate(alat::VectorOneVariableInterface* u, const solvers::FunctionInterface& function){_notWritten("interpolate");}
void FemInterface::toP1(alat::VectorOneVariableInterface* uc1, const alat::VectorOneVariableInterface* u){_notWritten("toP1");}
void FemInterface::fromP1(alat::VectorOneVariableInterface* u, const alat::VectorOneVariableInterface* uc1){_notWritten("fromP1");}
void FemInterface::computeErrors(int iK, solvers::ErrorsMap& errormaps, const arma::mat& uloc, const solvers::FunctionInterface& exactsolutions){_notWritten("computeErrors");}
bool FemInterface::noIntegration() const {return false;}

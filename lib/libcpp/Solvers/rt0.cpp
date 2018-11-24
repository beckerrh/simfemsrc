#include  "Alat/matrixonevariableinterface.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Solvers/rt0.hpp"
#include  "Solvers/integrationformulapoint.hpp"
#include  "Solvers/integrationformulaline.hpp"
#include  "Solvers/integrationformulatetrahedral.hpp"
#include  "Solvers/integrationformulatriangle.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
RT0::~RT0() {}
RT0::RT0(): solvers::Fem(){}
RT0::RT0( const RT0& c1): solvers::Fem(c1)
{
  (*this).operator=(c1);
}

RT0& RT0::operator=( const RT0& c1)
{
  solvers::Fem::operator=(c1);
  return *this;
}
std::string RT0::getClassName() const
{
  return "RT0";
}
std::unique_ptr<FemInterface> RT0::clone() const
{
  return std::unique_ptr<solvers::FemInterface>(new RT0(*this));
}
solverEnums::fem::femtype RT0::getType() const {return solverEnums::fem::RT0;}

/*--------------------------------------------------------------------------*/
int RT0::getN() const
{
  return _mesh->getNSides();
}
int RT0::getNPerCell(int iK) const
{
  return _mesh->getNSidesPerCell();
}
void RT0::indicesOfCell(int iK, alat::armaivec& indices) const
{
  int n = _mesh->getNSidesPerCell();
  assert(indices.size()==n);
  for(int ii=0;ii<n;ii++)
  {
    indices[ii] = _meshinfo->sides_of_cells(ii, iK);
  }
}
void RT0::setCell(int iK){_notWritten("setCell");}
void RT0::setCellBdry(int iK, int iS, int iil){_notWritten("setCellBdry");}
const FemData& RT0::referencePoint(const alat::Node& vhat, double weight){_notWritten("referencePoint");}
const FemData& RT0::referencePointBdry(const alat::Node& vhat, double weight){_notWritten("referencePointBdry");}

void RT0::strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const{_notWritten("strongDirichlet");}
void RT0::strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const{_notWritten("strongDirichletZero");}
void RT0::strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const{_notWritten("strongDirichlet");}
void RT0::setCellIsBdry(arma::uvec& cellisbdry){_notWritten("setCellIsBdry");}
const arma::uvec& RT0::getDofIsBdry() const{_notWritten("getNodeIsBdry");}
void RT0::setIsi(int iK){_notWritten("setIsi");}

/*--------------------------------------------------------------------------*/
void RT0::interpolate(alat::VectorOneVariableInterface* u, const solvers::FunctionInterface& function)
{
  alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
  int nloccell = getNPerCell();
  int nsides = _meshinfo->nsides;
  assert(uv->n()==nsides);
  assert(_ncomp==1);
  assert(uv->ncomp()==_ncomp);
  uv->fill(arma::fill::zeros);
  alat::armavec uinterpol(3,arma::fill::zeros);
  for(int iS=0; iS<nsides;iS++)
  {
    alat::Node xS = _mesh->getNodeOfSide(iS);
    function(uinterpol, xS.x(), xS.y(), xS.z(), _meshinfo->dim);
    // std::cerr << "uinterpol="<<uinterpol<< " normal="<<_meshinfo->normals.col(iS);
    (*uv)[iS] = arma::dot(_meshinfo->normals.col(iS), uinterpol);///arma::norm(_meshinfo->normals.col(iS));
  }
  // std::cerr << "uv="<<*uv;
}

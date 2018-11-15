#ifndef __Solvers_FemInterface_hpp
#define __Solvers_FemInterface_hpp

#include  "Alat/interfacebase.hpp"
#include  "Solvers/integrationformulainterface.hpp"
#include  "Solvers/applicationinterface.hpp"
#include  "Solvers/solverinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class Node;
  class MatrixOneVariableInterface;
  class MatrixAllVariables;
  class SparsityPatternSoft;
  class VectorOneVariableInterface;
}
namespace mesh
{
  class MeshUnitInterface;
}
namespace solvers
{
  class MeshInfo;
}
namespace solvers
{
  struct FemData
  {
    double J, x, y, z, weight, G;
    arma::vec phi, u, normal, uB, uI;
    arma::mat dphi, ugrad, uBgrad, uIgrad;
    arma::uvec isI;
    int iil, ncomp, nlocal;
    arma::mat mass, laplace;
    arma::vec mass_lumped;
  };
  class FemInterface : public virtual alat::InterfaceBase
  {
  protected:
    const mesh::MeshUnitInterface* _mesh;
    const MeshInfo* _meshinfo;
    int _ivar, _ncomp, _dim;
    virtual void initData() {}

  public:
    ~FemInterface();
    FemInterface();
    FemInterface( const FemInterface& feminterface);
    FemInterface& operator=( const FemInterface& feminterface);
    std::string getClassName() const;
    virtual std::unique_ptr<FemInterface> clone() const=0;
    virtual std::string getInfo() const;
    virtual solverEnums::fem::femtype getType() const=0;

    virtual void initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp);
    virtual int getN() const=0;
    virtual int getNPerCell(int iK=-1) const;
    virtual int getNcomp() const;
    virtual void indicesOfCell(int iK, alat::armaivec& indices) const;
    virtual const solvers::IntegrationFormulaInterface* getFormula() const;
    virtual const solvers::IntegrationFormulaInterface* getFormulaErrors() const;
    virtual const solvers::IntegrationFormulaInterface* getFormulaBdry() const;
    virtual const solvers::IntegrationFormulaInterface* getFormulaRhs() const;
    virtual void setCell(int iK);
    virtual void setCellBdry(int iK, int iS, int iil);
    virtual const FemData& referencePoint(const alat::Node& vhat, double weight);
    virtual const FemData& referencePointWithData(const alat::Node& vhat, double weight, const arma::mat& uloc);
    virtual const FemData& referencePointBdry(const alat::Node& vhat, double weight);
    virtual const FemData& referencePointBdryWithData(const alat::Node& vhat, double weight, const arma::mat& uloc);
    virtual const FemData& referencePointBdryCellWithData(const alat::Node& vhat, double weight, const arma::mat& uloc);
    virtual void setVectorIndices(int iK, alat::armaimat& vec_i)const;

    virtual void computeGrad(arma::mat& ugrad, const arma::mat& uloc) const;
    virtual void computeFunction(arma::vec& u, const arma::mat& uloc) const;
    virtual const solvers::FemData& getFemdata() const;
    virtual void setCellIsBdry(arma::uvec& cellisbdry);
    virtual void setIsi(int iK);

    virtual void strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const;
    virtual void strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const;
    virtual void strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const;
    virtual void interpolate(alat::VectorOneVariableInterface* u, const solvers::FunctionInterface& function);
    virtual bool canInterpolateToP1()const;
    virtual void toP1(alat::VectorOneVariableInterface* uc1, const alat::VectorOneVariableInterface* u);
    virtual void fromP1(alat::VectorOneVariableInterface* u, const alat::VectorOneVariableInterface* uc1);
    virtual void computeErrors(int iK, solvers::ErrorsMap& errormaps, const arma::mat& uloc, const solvers::FunctionInterface& exactsolutions);
    virtual bool noIntegration() const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

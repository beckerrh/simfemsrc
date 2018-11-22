#ifndef __Solvers_P2_hpp
#define __Solvers_P2_hpp

#include  "Solvers/fem.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  class P2 : public virtual solvers::Fem
  {
  protected:
    arma::vec _phip1;
    arma::mat _phip1grad;

    int _iil, _iK;
    alat::Node _vhat;
    arma::uvec _dofisbdry;
    void computeMatrices(int iK);
    std::unique_ptr<solvers::IntegrationFormulaInterface> newFormula();
    std::unique_ptr<solvers::IntegrationFormulaInterface> newFormulaErrors();
    std::unique_ptr<solvers::IntegrationFormulaInterface> newFormulaBdry();

  public:
    ~P2();
    P2();
    P2( const P2& p1);
    P2& operator=( const P2& p1);
    std::string getClassName() const;
    std::unique_ptr<solvers::FemInterface> clone() const;
    solverEnums::fem::femtype getType() const;
    void initData();

    int getN() const;
    int getNPerCell(int iK=-1) const;
    void indicesOfCell(int iK, alat::armaivec& indices) const;
    void setCell(int iK);
    void setCellBdry(int iK, int iS, int iil);
    const FemData& referencePoint(const alat::Node& vhat, double weight);
    const FemData& referencePointBdry(const alat::Node& vhat, double weight);

    void strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const;
    void strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const;
    void strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const;
    void setCellIsBdry(arma::uvec& cellisbdry);
    const arma::uvec& getDofIsBdry() const;
    void setIsi(int iK);
  };
}

/*--------------------------------------------------------------------------*/
#endif

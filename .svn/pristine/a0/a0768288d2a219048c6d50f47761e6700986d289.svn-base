#ifndef __Solvers_CR1_hpp
#define __Solvers_CR1_hpp

#include  "Solvers/fem.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  class CR1 : public virtual solvers::Fem
  {
  protected:
    int _iil;
    alat::Node _vhat;
    arma::uvec _dofisbdry;
    void computeMatrices(int iK);

  public:
    ~CR1();
    CR1();
    CR1( const CR1& p1);
    CR1& operator=( const CR1& p1);
    std::string getClassName() const;
    std::unique_ptr<solvers::FemInterface> clone() const;
    solverEnums::fem::femtype getType() const;

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

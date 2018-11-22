#ifndef __Solvers_RT0_hpp
#define __Solvers_RT0_hpp

#include  "Solvers/fem.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  class RT0 : public virtual solvers::Fem
  {
  protected:

  public:
    ~RT0();
    RT0();
    RT0( const RT0& p1);
    RT0& operator=( const RT0& p1);
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
    void interpolate(alat::VectorOneVariableInterface* u, const solvers::FunctionInterface& function) ;
    bool canInterpolateToP1()const{return false;}
  };
}

/*--------------------------------------------------------------------------*/
#endif

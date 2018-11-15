#ifndef __Solvers_Fem_hpp
#define __Solvers_Fem_hpp

#include  "Solvers/feminterface.hpp"
#include  "Solvers/meshunitwithdatainterface.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  class Fem : public virtual solvers::FemInterface
  {
  protected:
    std::shared_ptr<solvers::IntegrationFormulaInterface> _formula;
    std::shared_ptr<solvers::IntegrationFormulaInterface> _formulaerrors;
    std::shared_ptr<solvers::IntegrationFormulaInterface> _formulabdry;
    std::shared_ptr<solvers::IntegrationFormulaInterface> _formularhs;
    FemData _femdata;
    arma::vec _trafob;
    arma::mat _trafoA;
    virtual std::unique_ptr<solvers::IntegrationFormulaInterface> newFormula();
    virtual std::unique_ptr<solvers::IntegrationFormulaInterface> newFormulaErrors();
    virtual std::unique_ptr<solvers::IntegrationFormulaInterface> newFormulaBdry();
    virtual std::unique_ptr<solvers::IntegrationFormulaInterface> newFormulaRhs();
    void _computeData(const arma::mat& uloc);

  public:
    ~Fem();
    Fem();
    Fem(const Fem& fem);
    Fem& operator=( const Fem& fem);
    std::string getClassName() const;

    void initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp);
    const solvers::IntegrationFormulaInterface* getFormula() const;
    const solvers::IntegrationFormulaInterface* getFormulaErrors() const;
    const solvers::IntegrationFormulaInterface* getFormulaBdry() const;
    const solvers::IntegrationFormulaInterface* getFormulaRhs() const;
    const FemData& referencePointWithData(const alat::Node& vhat, double weight, const arma::mat& uloc);
    const FemData& referencePointBdryWithData(const alat::Node& vhat, double weight, const arma::mat& uloc);
    const FemData& referencePointBdryCellWithData(const alat::Node& vhat, double weight, const arma::mat& uloc);
    void computeGrad(arma::mat& ugrad, const arma::mat& uloc) const;
    void computeFunction(arma::vec& u, const arma::mat& uloc) const;
    const solvers::FemData& getFemdata() const;
    void setVectorIndices(int iK, alat::armaimat& vec_i)const;

    void toP1(alat::VectorOneVariableInterface* uc1, const alat::VectorOneVariableInterface* u);
    void fromP1(alat::VectorOneVariableInterface* u, const alat::VectorOneVariableInterface* uc1);
    void interpolate(alat::VectorOneVariableInterface* u, const solvers::FunctionInterface& function) ;
    void computeErrors(int iK, solvers::ErrorsMap& errormaps, const arma::mat& uloc, const solvers::FunctionInterface& exactsolutions);
  };
}

/*--------------------------------------------------------------------------*/
#endif

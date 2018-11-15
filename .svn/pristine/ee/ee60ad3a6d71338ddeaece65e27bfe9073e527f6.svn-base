#ifndef __Solvers_PdePartP1_hpp
#define __Solvers_PdePartP1_hpp

#include  "pdepartinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  class PdePartP1 : public PdePartInterface
  {
  protected:
    int _ivar, _ncomp;
    void _laplace(int iK, arma::mat& mat) const;
    void _massMatrix(int iK, arma::mat& mat, bool unscaled=false) const;

  public:
    ~PdePartP1();
    PdePartP1(alat::StringList vars);
    PdePartP1( const PdePartP1& pdepartc1);
    PdePartP1& operator=( const PdePartP1& pdepartc1);
    std::string getClassName() const;

    const solver_options::pdepart::opts setOptions()
    {
      return solver_options::pdepart::cell;
    }

    void setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters);

    void computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeRhsBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeResidualBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeMatrixCell(int iK, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const;
    void computeMatrixBdry(int color, int iK, int iS, int iil, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

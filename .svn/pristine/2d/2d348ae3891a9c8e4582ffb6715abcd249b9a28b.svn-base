#ifndef __Solvers_PdePartWithIntegration_hpp
#define __Solvers_PdePartWithIntegration_hpp

#include  "pdepartinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  class PdePartWithIntegration : public PdePartInterface
  {
  protected:

  public:
    ~PdePartWithIntegration();
    PdePartWithIntegration(alat::StringList vars);
    PdePartWithIntegration( const PdePartWithIntegration& pdepart);
    PdePartWithIntegration& operator=( const PdePartWithIntegration& pdepart);
    std::string getClassName() const;

    void rhsCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const;
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

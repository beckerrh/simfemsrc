#ifndef ___NewNitscheIntegration_hpp
#define ___NewNitscheIntegration_hpp

#include  "traditionalintegration.hpp"

/*--------------------------------------------------------------------------*/
class NewNitscheIntegration : public TraditionalIntegration
{
protected:
  mutable arma::mat _udirloc, _udirgrad;
  mutable arma::vec _udir;

public:
  ~NewNitscheIntegration();
  NewNitscheIntegration(alat::StringList vars);
  NewNitscheIntegration( const NewNitscheIntegration& nitscheintegration);
  NewNitscheIntegration& operator=( const NewNitscheIntegration& nitscheintegration);
  std::string getClassName() const;

  const solver_options::pdepart::opts setOptions();
  void rhsBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const;
  void residualBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const;
  void matrixBdry(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const;
  void prepareRhsCellBdry(int iK) const;
  void rhsBdryCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const;
  void residualBdryCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const;
  void matrixBdryCell(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const;
};

/*--------------------------------------------------------------------------*/
#endif

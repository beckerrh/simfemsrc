#ifndef ___NitscheIntegration_hpp
#define ___NitscheIntegration_hpp

#include  "traditionalintegration.hpp"

/*--------------------------------------------------------------------------*/
class NitscheIntegration : public TraditionalIntegration
{
protected:
  double _gammaNitsche(const solvers::FemData& fem)const;

public:
  ~NitscheIntegration();
  NitscheIntegration(alat::StringList vars);
  NitscheIntegration( const NitscheIntegration& nitscheintegration);
  NitscheIntegration& operator=( const NitscheIntegration& nitscheintegration);
  std::string getClassName() const;

  const solver_options::pdepart::opts setOptions();
  void rhsBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const;
  void residualBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const;
  void matrixBdry(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const;
};

/*--------------------------------------------------------------------------*/
#endif

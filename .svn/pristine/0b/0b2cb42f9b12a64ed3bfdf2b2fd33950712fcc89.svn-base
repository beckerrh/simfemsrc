#ifndef ___TraditionalIntegration_hpp
#define ___TraditionalIntegration_hpp

#include  "pdepart.hpp"
#include  "model.hpp"

/*--------------------------------------------------------------------------*/
class TraditionalIntegration : public PdePart
{
protected:

public:
  ~TraditionalIntegration();
  TraditionalIntegration(alat::StringList vars);
  TraditionalIntegration( const TraditionalIntegration& pdepartwithfemtraditional);
  TraditionalIntegration& operator=( const TraditionalIntegration& pdepartwithfemtraditional);
  std::string getClassName() const;

  const solver_options::pdepart::opts setOptions();
  void rhsCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const;
  void residualCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const;
  void matrixCell(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const;
};

/*--------------------------------------------------------------------------*/
#endif

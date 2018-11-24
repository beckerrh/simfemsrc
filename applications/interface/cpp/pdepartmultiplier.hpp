#ifndef ___PdePartMultiplier_hpp
#define ___PdePartMultiplier_hpp

#include  "pdepartcutinterface.hpp"
#include  "model.hpp"

/*--------------------------------------------------------------------------*/
class PdePartMultiplier : public PdePartCutInterface
{
protected:

public:
  ~PdePartMultiplier();
  PdePartMultiplier(alat::StringList vars);
  PdePartMultiplier( const PdePartMultiplier& pdepartwithfemtraditional);
  PdePartMultiplier& operator=( const PdePartMultiplier& pdepartwithfemtraditional);
  std::string getClassName() const;

  const solver_options::pdepart::opts setOptions();
  void additionCouplings(alat::Matrix<alat::SparsityPatternSoft>& sparsitypatternsoft)const;
  void computeMatrixGlobal(alat::MatrixAllVariables& A, const alat::VectorAllVariables& u)const;
  void computeResidualGlobal(alat::VectorAllVariables& r, const alat::VectorAllVariables& u)const;

  void computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const{}
  void computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const{}
  void computeMatrixCell(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const{}
};

/*--------------------------------------------------------------------------*/
#endif

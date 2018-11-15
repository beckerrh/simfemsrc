#ifndef ___PdePartIIM_hpp
#define ___PdePartIIM_hpp

#include  "pdepartcutinterface.hpp"
#include  "model.hpp"
#include  "p1iim.hpp"
#include  "cr1iim.hpp"

/*--------------------------------------------------------------------------*/
class PdePartIIM : public PdePartCutInterface
{
protected:
  double _gamma;
  std::shared_ptr<P1IIM> p1iim;
  std::shared_ptr<P1IIM> p1iimex;
  std::shared_ptr<CR1IIM> cr1iim;
  std::shared_ptr<CR1IIM> cr1iimex;
  void computeRhsCellCut(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
  void computeResidualCellCut(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
  void computeMatrixCellCut(int iK, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const;
  void computeMatrixIIM(int iK)const;
  void computeMatrixStabIIM(int iS, int iKin, int iKex)const;
  mutable arma::mat Linin, Linex, Lexin, Lexex, Last;

public:
  ~PdePartIIM();
  PdePartIIM(alat::StringList vars);
  PdePartIIM( const PdePartIIM& pdepartwithfemtraditional);
  PdePartIIM& operator=( const PdePartIIM& pdepartwithfemtraditional);
  std::string getClassName() const;

  const solver_options::pdepart::opts setOptions();
  void setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters);
  bool interiorsidecoupling(int iKin, int iKex) const;

  void computeResidualInteriorSide(int iS, int iKin, int iKex, solvers::PdePartData::vec& flocin, solvers::PdePartData::vec& flocex, const solvers::PdePartData::vec& ulocin, const solvers::PdePartData::vec& ulocex);
  void computeMatrixInteriorSide(int iS, int iKin, int iKex, solvers::PdePartData::mat& matinin, solvers::PdePartData::mat& matinex, solvers::PdePartData::mat& matexin, solvers::PdePartData::mat& matexex, solvers::PdePartData::imat& mat_i_in, solvers::PdePartData::imat& mat_j_in, solvers::PdePartData::imat& mat_i_ex, solvers::PdePartData::imat& mat_j_ex, const solvers::PdePartData::vec& ulocin, const solvers::PdePartData::vec& ulocex)const;
};

/*--------------------------------------------------------------------------*/
#endif

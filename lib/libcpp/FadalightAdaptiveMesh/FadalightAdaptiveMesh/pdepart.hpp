#ifndef ___PdePart_hpp
#define ___PdePart_hpp

#include  "pdepartcutinterface.hpp"
#include  "Solvers/rt0.hpp"
#include  "model.hpp"
#include  "p1cut.hpp"
#include  "p1hansbo.hpp"

/*--------------------------------------------------------------------------*/
class PdePart : public PdePartCutInterface
{
protected:
  double _gamma;
  std::shared_ptr<P1Cut> p1cut;
  std::shared_ptr<P1Hansbo> p1hansbo;
  std::string _method;
  mutable arma::mat Linin, Linex, Lexin, Lexex;

  void computeRhsCellCut(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
  void computeResidualCellCut(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
  void computeMatrixCellCut(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const;
  void computeProjectionNitsche(int iK)const;
  void computeProjectionNitscheStabC(int iK)const;
  void computeProjectionNitscheStabNC(int iK)const;
  void computeProjectionNewNitsche(int iK)const;
  void computeProjectionNewNitsche2(int iK)const;

public:
  ~PdePart();
  PdePart(alat::StringList vars);
  PdePart( const PdePart& pdepartwithfemtraditional);
  PdePart& operator=( const PdePart& pdepartwithfemtraditional);
  std::string getClassName() const;

  const solver_options::pdepart::opts setOptions();
  void setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters);
};

/*--------------------------------------------------------------------------*/
#endif

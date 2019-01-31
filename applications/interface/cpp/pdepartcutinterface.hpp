#ifndef ___PdePartCutInterface_hpp
#define ___PdePartCutInterface_hpp

#include  "Solvers/pdepartwithintegration.hpp"
#include  "model.hpp"

/*--------------------------------------------------------------------------*/
class PdePartCutInterface : public solvers::PdePartWithIntegration
{
protected:
  int _ivar, _ncomp, _nlocal;
  const Model* _localmodel;
  double _kin, _kex;

  const alat::armaivec* _cutcells, *_cutedges, *_celliscut, *_edgeiscut, *_cutnodes, *_nodeiscut, *_cutnodesisin;
  const alat::armavec* _cutcoeff;
  const alat::armaimat* _nodesofcutcellsisin, *_nodesofcutcells;
  const arma::mat* _measuresofcutcells, *_normalsofcutcells;
  const arma::mat* _cinofcutcells, *_cexofcutcells, *_cofinofcutcells, *_cofexofcutcells;

  alat::armaimat indicesofcutedge;
  alat::armamat coefsofcutedge;

  virtual void computeRhsCellCut(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const{}
  virtual void computeResidualCellCut(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const{}
  virtual void computeMatrixCellCut(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const{}

public:
  ~PdePartCutInterface();
  PdePartCutInterface(alat::StringList vars);
  PdePartCutInterface( const PdePartCutInterface& pdepartwithfemtraditional);
  PdePartCutInterface& operator=( const PdePartCutInterface& pdepartwithfemtraditional);
  std::string getClassName() const;

  const solver_options::pdepart::opts setOptions();
  void setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters);

  void computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
  void computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
  void computeMatrixCell(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const;
};

/*--------------------------------------------------------------------------*/
#endif

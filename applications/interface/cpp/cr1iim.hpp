#ifndef ___CR1IIM_hpp
#define ___CR1IIM_hpp

#include  "Solvers/p1.hpp"

/*--------------------------------------------------------------------------*/
class CR1IIM : public solvers::P1
{
protected:
  double _kin, _kex;
  arma::mat Iin, Iex, Bin, Bex, P, laplaceP1;
  alat::armavec phi_P1;
  arma::mat dphi_P1;

public:
  ~CR1IIM();
  CR1IIM();
  CR1IIM( const CR1IIM& P1cut);
  CR1IIM& operator=( const CR1IIM& P1cut);
  std::string getClassName() const;
  const alat::armaivec* _cutcells, *_cutedges, *_celliscut, *_edgeiscut, *_cutnodes;
  const alat::armavec* _cutcoeff;
  const alat::armaimat* _nodesofcutcellsisin, *_nodesofcutcells;
  const arma::mat* _measuresofcutcells, *_normalsofcutcells;
  const arma::mat* _cinofcutcells, *_cexofcutcells, *_cofinofcutcells, *_cofexofcutcells;
  int iKcut;
  bool iKisin;
  arma::mat femcoefin, femcoefex, laplaceCut;

  int getN() const;
  int getNPerCell(int iK=-1) const;
  void indicesOfCell(int iK, alat::armaivec& indices) const;
  void setCoefs(double kin, double kex);
  void setCell(int iK);
  const solvers::FemData& referencePoint(const alat::Node& vhat, double weight);

  void initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp);

  void computeCutFem(int iK);

  void initCutInterface(const mesh::MeshUnitInterface* mesh);
  void computeErrors(int iK, solvers::ErrorsMap& errormaps, const alat::armavec& uloc, const solvers::FunctionInterface& exactsolutions);

  void strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const;
  void strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const;
  void strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const;
};

/*--------------------------------------------------------------------------*/
#endif

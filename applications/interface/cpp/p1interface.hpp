#ifndef ___P1Interface_hpp
#define ___P1Interface_hpp

#include  "Solvers/feminterface.hpp"

/*--------------------------------------------------------------------------*/
class P1Interface : public solvers::FemInterface
{
protected:
  alat::armaivec getBdryNodes(int i, const alat::armaimat& cells_on_bdry)const;
  double _kin, _kex;

public:
  ~P1Interface();
  P1Interface();
  P1Interface( const P1Interface& P1cut);
  P1Interface& operator=( const P1Interface& P1cut);
  std::string getClassName() const;
  const alat::armaivec* _cutcells, *_cutedges, *_celliscut, *_edgeiscut, *_cutnodes, *_nodeiscut, *_cutnodesisin;
  const arma::vec* _cutcoeff;
  const alat::armaimat* _nodesofcutcellsisin, *_nodesofcutcells;
  const arma::mat* _measuresofcutcells, *_normalsofcutcells;
  const arma::mat* _cinofcutcells, *_cexofcutcells, *_cofinofcutcells, *_cofexofcutcells;

  alat::armaivec nodeofcutedgelim, nodeofcutedgeto;
  alat::armaimat indicesofcutedge;
  alat::armamat coefsofcutedge, coefsofcutedgeweights;

  void initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp);
  bool noIntegration() const;

  std::unique_ptr<solvers::FemInterface> clone() const;
  solverEnums::fem::femtype getType() const;

  int getN() const;
  int getNPerCell(int iK) const;
  void indicesOfCell(int iK, alat::armaivec& indices) const;
  void setCell(int iK);

  void initCutInterface(const mesh::MeshUnitInterface* mesh);
  bool canInterpolateToP1()const{return false;}
  void strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const {}
  void strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const {}
  void strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const {}

};

/*--------------------------------------------------------------------------*/
#endif

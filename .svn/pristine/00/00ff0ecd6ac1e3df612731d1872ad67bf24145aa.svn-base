#ifndef ___P1Hansbo_hpp
#define ___P1Hansbo_hpp

#include  "Solvers/p1.hpp"

/*--------------------------------------------------------------------------*/
class P1Hansbo : public solvers::P1
{
protected:
  alat::armaivec getBdryNodes(int i, const alat::armaimat& cells_on_bdry)const;
  double _kin, _kex;

public:
  ~P1Hansbo();
  P1Hansbo();
  P1Hansbo( const P1Hansbo& P1cut);
  P1Hansbo& operator=( const P1Hansbo& P1cut);
  std::string getClassName() const;
  const alat::armaivec* _cutcells, *_cutedges, *_celliscut, *_edgeiscut, *_cutnodes, *_nodeiscut, *_cutnodesisin;
  const arma::vec* _cutcoeff;
  const alat::armaimat* _nodesofcutcellsisin, *_nodesofcutcells;
  const arma::mat* _measuresofcutcells, *_normalsofcutcells;
  const arma::mat* _cinofcutcells, *_cexofcutcells, *_cofinofcutcells, *_cofexofcutcells;

  alat::armaivec indin, index;

  void initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp);
  void setCoefs(double kin, double kex);

  int getN() const;
  int getNPerCell(int iK) const;
  void indicesOfCell(int iK, alat::armaivec& indices) const;
  void setCell(int iK);

  void initCutInterface(const mesh::MeshUnitInterface* mesh);
  void toP1(alat::VectorOneVariableInterface* uP1, const alat::VectorOneVariableInterface* u);
  void fromP1(alat::VectorOneVariableInterface* u, const alat::VectorOneVariableInterface* uP1);
  void computeErrors(int iK, solvers::ErrorsMap& errormaps, const arma::mat& uloc, const solvers::FunctionInterface& exactsolutions);
  void strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const;
  void strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const;
  void strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const;

};

/*--------------------------------------------------------------------------*/
#endif

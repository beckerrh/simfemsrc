#ifndef ___P1IIM_hpp
#define ___P1IIM_hpp

#include  "Solvers/p1.hpp"

/*--------------------------------------------------------------------------*/
class P1IIM : public solvers::P1
{
protected:
  double _kin, _kex;
  arma::mat Iin, Iex, Bin, Bex, C, P, laplaceP1;

public:
  ~P1IIM();
  P1IIM();
  P1IIM( const P1IIM& P1cut);
  P1IIM& operator=( const P1IIM& P1cut);
  std::string getClassName() const;
  const alat::armaivec* _cutcells, *_cutedges, *_celliscut, *_edgeiscut, *_cutnodes;
  const alat::armavec* _cutcoeff;
  const alat::armaimat* _nodesofcutcellsisin, *_nodesofcutcells;
  const arma::mat* _measuresofcutcells, *_normalsofcutcells;
  const arma::mat* _cinofcutcells, *_cexofcutcells, *_cofinofcutcells, *_cofexofcutcells;
  int iKcut;
  bool iKisin;
  arma::mat femcoefin, femcoefex, laplaceCut;

  void setCoefs(double kin, double kex);

  void initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp);

  void setCell(int iK);

  void computeCutFem(int iK);

  void initCutInterface(const mesh::MeshUnitInterface* mesh);
  void computeErrors(int iK, solvers::ErrorsMap& errormaps, const alat::armavec& uloc, const solvers::FunctionInterface& exactsolutions);

};

/*--------------------------------------------------------------------------*/
#endif

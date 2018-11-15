#ifndef __Mesh_CutInterface_h
#define __Mesh_CutInterface_h

#include  "Alat/vector.hpp"
#include  "Mesh/geometryobject.hpp"

/*--------------------------------------------------------------------------*/

namespace alat
{
  class VectorOneVariable;
}

namespace mesh
{
  class CutInterfaceConstructor : public GeometryConstructorInterface{
  public: virtual const alat::VectorOneVariable* getPhi() const=0;
  std::string getClassName() const {return "CutInterfaceConstructor";}
};

  class CutInterface : public GeometryObject
  {
protected:
    alat::armaivec _cutcells, _cutedges, _cutnodes, _celliscut, _edgeiscut, _nodeiscut, _cutnodesisin;
    arma::vec _cutcoeff;
    alat::armaimat _nodesofcutcellsisin, _nodesofcutcell;
    arma::mat _measuresofcutcells, _normalsofcutcells;
    arma::mat _cinofcutcells, _cexofcutcells, _cofinofcutcells, _cofexofcutcells;

    void _computeMeasuresOfCutCell(int iK, const mesh::MeshUnitInterface* mesh, arma::subview_col<double> measuresofcutcell, arma::subview_col<double> normalsofcutcells, arma::subview_col<double> cinofcutcells, arma::subview_col<double> cexofcutcells, const arma::mat& innodes, const arma::mat& exnodes, const arma::mat& addnodes, double moc, const arma::ivec& innodesind, const arma::ivec& exnodesind, const arma::ivec& addnodesedges);
    double _surface(arma::subview_col<double> u, arma::subview_col<double> v, arma::subview_col<double> w)const;
    double _surface(const arma::mat& addnodes)const;
    void _computeCutTet(double& volin, double& volex, arma::subview_col<double> xcin, arma::subview_col<double> xcex, const arma::vec& normal, int iK, const mesh::MeshUnitInterface* mesh, const arma::mat& innodes, const arma::mat& exnodes, const arma::mat& addnodes, const arma::ivec& innodesind, const arma::ivec& exnodesind, const arma::ivec& addnodesedges)const;

public:
    ~CutInterface();
    CutInterface();
    CutInterface( const CutInterface& geometryobject);
    CutInterface& operator=( const CutInterface& geometryobject);
    std::string getClassName() const;
    std::unique_ptr<GeometryObject> clone() const;

    const alat::armaivec& getCutCells() const;
    const alat::armaivec& getCutEdges() const;
    const alat::armaivec& getCellIsCut() const;
    const alat::armaivec& getEdgeIsCut() const;
    const alat::armaivec& getNodeIsCut() const;
    const alat::armaivec& getCutNodeIsIn() const;
    const alat::armaivec& getCutNodes() const;
    const arma::vec& getCutCoeff() const;
    const alat::armaimat& getNodesOfCutCellsIsIn() const;
    const alat::armaimat& getNodesOfCutCells() const;
    const arma::mat& getMeasureOfCutCells() const;
    const arma::mat& getNormalsOfCutCells() const;
    const arma::mat& getCInOfCutCells() const;
    const arma::mat& getCExOfCutCells() const;
    const arma::mat& getCofInOfCutCells() const;
    const arma::mat& getCofExOfCutCells() const;

    alat::armaivec getSizes() const;
    void setSizes(alat::armaivec::const_iterator sizes);
    void send(int neighbor, int tag) const;
    void recv(int neighbor, int tag);
    void loadH5(const arma::hdf5_name& spec);
    void saveH5(const arma::hdf5_name& spec) const;
    void construct(const mesh::MeshUnitInterface* mesh, const GeometryConstructorInterface* geometryconstructor);
  };
}

/*--------------------------------------------------------------------------*/

#endif

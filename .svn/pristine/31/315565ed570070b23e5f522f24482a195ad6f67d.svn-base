
#ifndef __Mesh_MeshUnit_h
#define __Mesh_MeshUnit_h

#include  "Alat/map.hpp"
#include  "Mesh/geometryobjectsconstructor.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Mesh/sidesandcells.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Mesh/interfacemeshinfo.hpp"
#include  "Mesh/meshvisitorinterface.hpp"
#include  "Alat/armadillo.hpp"

/*---------------------------------------------------------*/
namespace alat
{
  class SparsityPattern;
}
namespace mesh
{
  class MeshUnit : public virtual mesh::MeshUnitInterface
  {
  public:
    bool _init_called;
    int n_dim, n_nodes_per_cell, n_sides_per_cell, n_nodes_per_side, n_edges_per_cell;

  protected:
    std::shared_ptr<const mesh::MeshVisitorInterface> _visitor;
    NodesAndNodesOfCells _nodesandnodesofcells;
    SidesAndCells _sidesandcells;
    EdgesAndCells _edgesandcells;
    BoundaryInformationMap _boundaryinformationmap;
    GeometryObjectsConstructor _geometryobjects_constructor;
    GeometryObjects _geometryobjects;

  public:
    ~MeshUnit();
    MeshUnit();
    MeshUnit(const MeshUnit& meshunit);
    MeshUnit& operator=(const MeshUnit& meshbase);
    void init(std::shared_ptr<const mesh::MeshVisitorInterface> visitor);

    const mesh::MeshVisitorInterface* getVisitor() const;
    std::string getClassName() const;
    std::string getInfo() const;
    std::string getXdmfTopoType() const;
    int getVtkCellType() const;
    int getVtkSideType() const;
    int getGmshCellType() const;
    int getGmshSideType() const;
    int getGmshSideSideType() const;
    int getGmshSideSideSideType() const;
    alat::armaivec _getSideOfCell(int i, int ii) const;
    alat::armaivec _getEdgeOfCell(int i, int ii) const;

    const BoundaryInformationMap& getBoundaryInformationMap()const;
    BoundaryInformationMap& getBoundaryInformationMap();
    const BoundaryInformation& getBoundaryInformation(int color) const;
    BoundaryInformation& getBoundaryInformation(int color);
    const NodesAndNodesOfCells& getNodesAndNodesOfCells() const;
    NodesAndNodesOfCells& getNodesAndNodesOfCells();
    const SidesAndCells& getSidesAndCells() const;
    SidesAndCells& getSidesAndCells();
    const EdgesAndCells& getEdgesAndCells() const;
    EdgesAndCells& getEdgesAndCells();

    const GeometryObjects& getGeometryObjects() const;
    GeometryObjects& getGeometryObjects();
    bool geometryObjectExists(meshEnums::geomobjtype type) const;
    std::shared_ptr<const GeometryObject> getGeometryObject(meshEnums::geomobjtype type) const;
    std::shared_ptr<GeometryObject> getGeometryObject(meshEnums::geomobjtype type);
    void addGeometryObject(meshEnums::geomobjtype type, const GeometryConstructorInterface* geometryconstructor=&trivialgeometryconstructor);

    int getDimension() const;
    int getNNodesPerCell() const;
    int getNNodesPerSide() const;
    int getNEdgesPerCell() const;
    int getNSidesPerCell() const;
    int getNNodes() const;
    int getNCells() const;
    int getNSides() const;
    int getNEdges() const;

    const arma::mat&  getNodes() const;
    arma::mat& getNodes();
    const alat::armaimat& getCells() const;
    alat::armaimat& getCells();
    const alat::armaimat& getSides() const;
    alat::armaimat& getSides();
    const alat::armaimat& getSidesOfCells() const;
    alat::armaimat& getSidesOfCells();
    const alat::armaimat& getCellsOfSides() const;
    alat::armaimat& getCellsOfSides();

    const MeshUnitInterface::Cell getCell(int i) const;
    const MeshUnitInterface::Side getSide(int i) const;
    alat::Node getNodeOfCell(int iK) const;
    alat::Node getNodeOfSide(int iS) const;

    void saveH5(const arma::hdf5_name& spec) const;
    void loadH5(const arma::hdf5_name& spec);
    void writeVtk(std::string filename) const;
    void writeBoundaryVtk(std::string filename) const;

    alat::SparsityPattern getSizes() const;
    void setSizes(const alat::SparsityPattern& sparsitypattern);
    void sendRecv(std::shared_ptr<mesh::MeshUnitInterface> meshunit, int neighbor) const;
  };
}

/*---------------------------------------------------------*/

#endif

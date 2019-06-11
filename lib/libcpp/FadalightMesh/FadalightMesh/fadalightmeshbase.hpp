#ifndef __FadalightMesh_FadalightMeshBase_h
#define __FadalightMesh_FadalightMeshBase_h

#include  "FadalightMesh/boundaryinfo.hpp"
#include  "FadalightMesh/curvedinteriorsideinfo.hpp"
#include  "FadalightMesh/geometryobjectsconstructorinterface.hpp"
#include  "FadalightMesh/meshinterface.hpp"
#include  "curvedboundaryinformation.hpp"
#include  <sstream>

/*---------------------------------------------------------*/

namespace FadalightMesh
{
  template<int N>
  class CellBase : public alat::FixArray<N, int>
  {
public:
    int node(int i) const;
  };

  /*---------------------------------------------------------*/


  template<int DIM, int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
  class FadalightMeshBase : public virtual FadalightMesh::MeshInterface
  {
public:
    typedef CellBase<NODESPERCELL> Cell;
    typedef CellBase<NODESPERSIDE> Side;
    typedef CellBase<SIDESPERCELL> SideCell;
    typedef CellBase<2> CellSide;
    typedef alat::Map<Side, int> BoundarySideToColor;

protected:
    mutable std::string _visutype;
    static Side _dummyside;
    std::string _infilename;
    alat::Vector<alat::Node>   _nodes;
    alat::Vector<Cell>   _cells;
    alat::Vector<Side>   _sides;
    alat::Vector<SideCell>   _sides_of_cells;
    alat::Vector<CellSide>   _cells_of_sides;
    // alat::SparsityPatternFixArray<2> _neighbors_of_cells;

    FadalightMesh::BoundaryInfo _boundaryinfo;
    FadalightMesh::CurvedInteriorSideInfo _curvedinteriorsideinfo;

    FadalightMesh::GeometryObjectsConstructorInterface* _geometryobjects_constructor;
    std::map<std::string, FadalightMesh::GeometryObject*>   _geometryobjects;

    //! this function is need to construct the side information from _cells
    virtual const Side& _getSideOfCell(int i, int ii) const;

    void checkGeoFile(std::string filename) const;

public:
    FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>( );
    FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>( const FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>&fadalightmeshbase );
    ~FadalightMeshBase<DIM, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>( );

    const std::string& getInfileName() const;
    void setVisuType(const std::string& visutype) const;
    int getNLevels() const;
    void setResolution(int level) const;
    int getResolution() const;
    void getNodesOfCell(int iK, alat::Vector<alat::Node>& F) const;

    // void reInit();
    void createGeometryObject(std::string name);
    bool geometryObjectExists(std::string name) const;
    const FadalightMesh::GeometryObject* getGeometryObject(std::string name) const;
    FadalightMesh::GeometryObject* getGeometryObject(std::string name);

    // fonctions d'accès au geometryobjects_constructor
    // const FadalightMesh::GeometryObjectsConstructorInterface* getGeometryObjectsConstructor() const;
    // FadalightMesh::GeometryObjectsConstructorInterface*& getGeometryObjectsConstructor();
    // void addGeometryObject(const std::string& name, FadalightMesh::GeometryObject* geo );
    // const std::map<std::string, FadalightMesh::GeometryObject*>&  getGeometryObjects() const;
    // std::map<std::string, FadalightMesh::GeometryObject*>&  getGeometryObjects();
    alat::Vector<alat::Node>&   getAllNodes();
    alat::Vector<Cell>&   getCells();
    alat::Vector<Side>&   getSides();
    alat::Vector<SideCell>&  getSidesOfCells();
    alat::Vector<CellSide>&  getCellsOfSides();
    // CurvedBoundaryInformation* getCurvedBoundaryInformation();
    // const CurvedBoundaryInformation* getCurvedBoundaryInformation() const;
    const alat::Vector<alat::Node>&   getNodes() const;
    const alat::Vector<Cell>&   getCells() const;
    const alat::Vector<Side>&   getSides() const;
    const alat::Vector<SideCell>&  getSidesOfCells() const;
    const alat::Vector<CellSide>&  getCellsOfSides() const;
    const Cell&   getCell(int i) const;
    const Side&   getSide(int i) const;
    const SideCell&  getSidesOfCell(int i) const;
    std::string getClassName() const;
    std::string getInfo() const;
    int getDimension() const;
    int getNNodesPerCell(int i) const;
    int getNNodesPerSide(int i) const;
    int getNSidesPerCell(int i) const;
    int getNNodes() const;
    int getNCells() const;
    int getNSides() const;
    int getNEdges() const;
    int getNodeIdOfCell(int iK, int ii) const;
    int getNodeIdOfSide(int iS, int ii) const;
    int getSideIdOfCell(int iK, int ii) const;
    int getCellIdOfSide(int iS, int ii) const;
    const alat::Node& getNode(int i) const;
    const alat::Node& getNodeOfSide(int is, int ii) const;
    const alat::Node& getNodeOfCell(int iK, int ii) const;
    alat::Node getNodeOfCell(int iK) const;
    alat::Node getNodeOfSide(int iS) const;

    int getLocalIndexOfSideInCell(int iK, int iS) const;
    const FadalightMesh::BoundaryInfo* getBoundaryInfo() const;
    FadalightMesh::BoundaryInfo* getBoundaryInfo();
    const FadalightMesh::CurvedInteriorSideInfo* getCurvedInteriorSideInfo() const;
    // FadalightMesh::CurvedInteriorSideInfo* getCurvedInteriorSideInfo();
    const alat::armaivec& getBoundaryColors() const;
    const alat::armaivec& getBoundarySides(int color) const;
    // bool cellIsCurved(int iK) const;
    int findNeightborHangingCells(int iK, int iS, alat::Node pt);
    void readFadalightMesh(const std::string& basefilename);
    void writeFadalightMesh(const std::string& basefilename, arma::file_type datatype = arma::arma_binary) const;
    // void writeEnsightGeometryObjectDescFile(const std::string& basefilename);
    void writeSimpleMesh(std::string filename, arma::file_type datatype = arma::arma_binary) const;

    ///! supposes that _cells is filled; fills _sides and _sides_of_cells.
    void constructSidesFromCells(BoundarySideToColor& bsides, int color_default = 0);
    void constructSidesFromCells(BoundarySideToColor& bsides, const BoundarySideToColor& icsides, int color_default = 0);
    // void computeCellNeighbours();
    // void computeCellConnectivity(alat::SparsityPattern& SPC) const;

    void writeVtk(std::string filename) const;
    // void writeEnsightGeometry(std::string filename) const;
    void writeMeshInfo(std::string filename) const;
    void writeH5(H5::H5File& file) const;

    void writeBoundaryVtk(std::string filename) const;
    // void addGeometryObjects(const std::string& filename, const alat::StringVector& names_of_geometry_objects, std::string datatype);
    bool isMultilevel() const;

    int localIndexOfNodeInCell(int iK, int in) const;
    int localIndexOfSideInCell(int iK, int is) const;
    int localIndexOfEdgeInCell(int iK, int ie) const;

    // const alat::SparsityPatternFixArray<2>& getCellNeighbours() const;
  };
}

/*---------------------------------------------------------*/

#endif

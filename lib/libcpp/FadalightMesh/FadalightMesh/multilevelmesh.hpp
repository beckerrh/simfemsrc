#ifndef __FadalightMesh_MultiLevelMesh_h
#define __FadalightMesh_MultiLevelMesh_h

#include  "FadalightMesh/meshinterface.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class MultiLevelMesh : public FadalightMesh::MeshInterface
  {
protected:
    mutable int _activelevel;
    std::string _type;
    alat::Vector<FadalightMesh::MeshInterface*> _meshes;

public:
    ~MultiLevelMesh();
    MultiLevelMesh(std::string type);
    MultiLevelMesh( const MultiLevelMesh& multilevelmesh);
    MultiLevelMesh& operator=( const MultiLevelMesh& multilevelmesh);
    std::string getClassName() const;
    MultiLevelMesh* clone() const;

    alat::Vector<FadalightMesh::MeshInterface*>& getMeshes();
    std::string getCellType() const;
    int getNLevels() const;
    const FadalightMesh::MeshInterface* getMesh(int level) const;
    FadalightMesh::MeshInterface* getMesh(int level);
    FadalightMesh::MeshInterface*& getMeshPointer(int level);
    std::string getInfo() const;
    int getDimension() const;
    void readFadalightMesh(const std::string& basefilename);
    void writeFadalightMesh(const std::string& basefilename, arma::file_type datatype = arma::arma_binary) const;
    const FadalightMesh::BoundaryInfo* getBoundaryInfo() const;
    void setResolution(int level) const;
    int getResolution() const;
    int getNNodes() const;
    int getNEdges() const;
    int getNSides() const;
    int getNCells() const;
    int getNNodesPerCell(int iK) const;
    int getNNodesPerSide(int iS) const;
    int getNSidesPerCell(int iK) const;
    int getNEdgesPerCell(int iK) const;
    int getNEdgesPerSide(int iS) const;
    int getNodeIdOfCell(int i, int ii) const;
    int getNodeIdOfSide(int i, int ii) const;
    int getNodeIdOfEdge(int i, int ii) const;
    int getSideIdOfCell(int i, int ii) const;
    int getEdgeIdOfCell(int i, int ii) const;
    int getEdgeIdOfSide(int i, int ii) const;
    int getCellIdOfSide(int i, int ii) const;
    int getLocalIndexOfSideInCell(int iK, int iS) const;
    void getNodesOfCell(int iK, alat::Vector<alat::Node>& F) const;
    void writeH5(H5::H5File& file) const;
    void writeMeshInfo(std::string filename) const;
    bool geometryObjectExists(std::string name) const;
    const FadalightMesh::GeometryObject* getGeometryObject(std::string name) const;
    bool isMultilevel() const;
    int getCouplingOffset(int iS) const;
    void setVisuType(const std::string& visutype) const;
 
    alat::Node getNodeOfCell(int iK) const;
    alat::Node getNodeOfSide(int iS) const;

    int localIndexOfNodeInCell(int iK, int in) const;
    int localIndexOfSideInCell(int iK, int is) const;
    int localIndexOfEdgeInCell(int iK, int ie) const;
    const alat::Node& getNodeOfSide(int is, int ii) const;
    const alat::Node& getNode(int i) const;
    FadalightMeshEnums::meshtype getType() const;
    // bool cellIsCurved(int iK) const;
    void getLocalIndicesOfSidesInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const;
    void getLocalIndicesOfSidesAndDiagonalsInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

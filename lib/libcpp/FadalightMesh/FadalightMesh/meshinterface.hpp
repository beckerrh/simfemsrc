#ifndef __FadalightMesh_MeshInterface_h
#define __FadalightMesh_MeshInterface_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/node.hpp"
#include  "Alat/sparsitypatternfixarray.hpp"
#include  "enums.hpp"

/*--------------------------------------------------------------------------*/
namespace H5
{
  class H5File;
}
namespace alat
{
  class StringVector;
}
namespace FadalightMesh
{
  class BoundaryInfo;
  class CurvedBoundaryInformation;
  class CurvedInteriorSideInfo;
  class GeometryObject;
  class GeometryObjectsConstructorInterface;
  class SparsityPattern;

  class MeshInterface : public alat::InterfaceBase
  {
private:
    static alat::Node _dummynode;
    GeometryObjectsConstructorInterface* dummy_geometryobjects_constructor;

protected:
    std::string getInterfaceName() const;

public:
    ~MeshInterface();
    MeshInterface();
    MeshInterface( const MeshInterface& meshinterface);
    MeshInterface& operator=( const MeshInterface& meshinterface);

    virtual const std::string& getInfileName() const;
    virtual void setResolution(int level) const;
    virtual int getResolution() const;
    virtual int getNLevels() const;
    virtual std::string getInfo() const = 0;
    virtual void getNodesOfCell(int iK, alat::Vector<alat::Node>& F) const;
    virtual int getDimension() const = 0;
    virtual alat::Node getNodeOfSide(int iS) const;
    virtual int getNodeIdOfSideOfCell(int iK, int iis, int ii) const;
    virtual alat::Node getNodeOfCell(int iK) const;
    virtual alat::Node getNodeOfEdge(int iK) const;
    virtual int getNNodes() const;
    virtual int getNEdges() const;
    virtual int getNSides() const;
    virtual int getNCells() const;
    virtual const alat::Node& getNode(int i) const;
    virtual const alat::Node& getNodeOfCell(int i, int ii) const;
    virtual const alat::Node& getNodeOfEdge(int i, int ii) const;
    virtual const alat::Node& getNodeOfSide(int i, int ii) const;

    virtual alat::Vector<alat::Node>& getAllNodes();

    virtual std::string getCellType() const;
    virtual int getNNodesPerCell(int iK) const;
    virtual int getNNodesPerSide(int iS) const;
    virtual int getNSidesPerCell(int iK) const;
    virtual int getNEdgesPerCell(int iK) const;
    virtual int getNEdgesPerSide(int iS) const;
    virtual const FadalightMesh::BoundaryInfo* getBoundaryInfo() const;
    virtual const FadalightMesh::CurvedInteriorSideInfo* getCurvedInteriorSideInfo() const ;
    virtual int getNodeIdOfCell(int i, int ii) const;
    virtual int getNodeIdOfSide(int i, int ii) const;
    virtual int getNodeIdOfEdge(int i, int ii) const;
    virtual int getSideIdOfCell(int i, int ii) const;
    virtual int getEdgeIdOfCell(int i, int ii) const;
    virtual int getEdgeIdOfSide(int i, int ii) const;
    virtual int getCellIdOfSide(int i, int ii) const;
    virtual int getLocalIndexOfSideInCell(int iK, int iS) const;
    virtual bool geometryObjectExists(std::string name) const;
    virtual const FadalightMesh::GeometryObject* getGeometryObject(std::string name) const;
    virtual FadalightMesh::GeometryObject* getGeometryObject(std::string name);
    virtual void createGeometryObject(std::string name);
    virtual const FadalightMesh::CurvedBoundaryInformation* getCurvedBoundaryInformation() const;
    virtual void read(std::string filename);
    virtual void write(std::string filename) const;
    virtual void readFadalightMesh(const std::string& basefilename);
    virtual void writeFadalightMesh(const std::string& basefilename, arma::file_type datatype = arma::arma_binary) const;
    virtual void writeSimpleMesh(std::string basefilename, arma::file_type datatype = arma::arma_binary) const;
    virtual void writeVtk(std::string filename) const;
    virtual void writeBoundaryVtk(std::string filename) const;
    virtual void writeMeshInfo(std::string filename) const;
    virtual void writeH5(H5::H5File& file) const;
    virtual void setVisuType(const std::string& visutype) const = 0;

    virtual int getVtkType() const;
    virtual int getBoundaryVtkType() const;
    virtual bool cellIsCurved(int iK) const;
    virtual int findNeightborHangingCells(int iK, int iS, alat::Node pt);
    virtual int getLocalNodeIndiceOfSide(int ii, int isl) const;
    virtual int getCouplingOffset(int iS) const;
    virtual bool isMultilevel() const = 0;

    virtual int localIndexOfNodeInCell(int iK, int in) const;
    virtual int localIndexOfSideInCell(int iK, int is) const;
    virtual int localIndexOfEdgeInCell(int iK, int ie) const;

    virtual FadalightMeshEnums::meshtype getType() const;
    virtual void getLocalIndicesOfSidesInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const;
    virtual void getLocalIndicesOfSidesAndDiagonalsInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

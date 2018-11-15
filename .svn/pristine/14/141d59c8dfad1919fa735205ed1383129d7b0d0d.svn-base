#ifndef __Mesh_MeshUnitInterface_h
#define __Mesh_MeshUnitInterface_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/map.hpp"
#include  "Alat/node.hpp"
#include  "Mesh/enums.hpp"
#include  "Mesh/boundaryinformation.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class GeometryObject;
  class MeshVisitorInterface;
  class NodesAndNodesOfCells;
  class EdgesAndCells;
  class SidesAndCells;

  class MeshUnitInterface : public virtual alat::InterfaceBase
  {
  public:
    typedef arma::subview_col<double> Normal;
//  #ifdef ARMA_64BIT_WORD
//     typedef arma::subview_col<arma::sword> Cell;
//     typedef arma::subview_col<arma::sword> Side;
// #else
    typedef arma::subview_col<int> Cell;
    typedef arma::subview_col<int> Side;
// #endif
    typedef alat::Map<int, BoundaryInformation> BoundaryInformationMap;
    typedef alat::Map<int,alat::armaimat>  BoundarySideToColorForSidesConstructor;
    typedef alat::Map<meshEnums::geomobjtype, std::shared_ptr<mesh::GeometryObject> >  GeometryObjects;

  protected:
    std::string getInterfaceName() const;

  public:
    ~MeshUnitInterface();
    MeshUnitInterface();
    MeshUnitInterface( const MeshUnitInterface& meshinterface);
    MeshUnitInterface& operator=( const MeshUnitInterface& meshinterface);

    virtual const mesh::MeshVisitorInterface* getVisitor() const=0;
    virtual std::string getInfo() const = 0;
    virtual std::string getXdmfTopoType() const=0;
    virtual int getVtkCellType() const=0;
    virtual int getVtkSideType() const=0;
    virtual int getGmshCellType() const=0;
    virtual int getGmshSideType() const=0;
    virtual int getGmshSideSideType() const=0;
    virtual int getGmshSideSideSideType() const=0;

    //! this function is used to construct the side information from _cells
    virtual alat::armaivec _getSideOfCell(int i, int ii) const;
    virtual alat::armaivec _getEdgeOfCell(int i, int ii) const;

    // BoundaryInformation
    virtual const BoundaryInformationMap& getBoundaryInformationMap()const=0;
    virtual BoundaryInformationMap& getBoundaryInformationMap()=0;
    virtual const BoundaryInformation& getBoundaryInformation(int color) const=0;
    virtual BoundaryInformation& getBoundaryInformation(int color)=0;
    virtual const NodesAndNodesOfCells& getNodesAndNodesOfCells() const=0;
    virtual NodesAndNodesOfCells& getNodesAndNodesOfCells()=0;
    virtual const SidesAndCells& getSidesAndCells() const=0;
    virtual SidesAndCells& getSidesAndCells()=0;
    virtual const EdgesAndCells& getEdgesAndCells() const=0;
    virtual EdgesAndCells& getEdgesAndCells()=0;
   // GeometryObjects
    virtual const GeometryObjects& getGeometryObjects() const=0;
    virtual GeometryObjects& getGeometryObjects()=0;
    virtual bool geometryObjectExists(meshEnums::geomobjtype type) const=0;
    virtual std::shared_ptr<const GeometryObject> getGeometryObject(meshEnums::geomobjtype type) const=0;
    virtual std::shared_ptr<GeometryObject> getGeometryObject(meshEnums::geomobjtype type)=0;
    virtual void addGeometryObject(meshEnums::geomobjtype type, const GeometryConstructorInterface* geometryconstructor=&trivialgeometryconstructor)=0;

    virtual int getDimension() const = 0;
    virtual int getNNodesPerCell() const=0;
    virtual int getNNodesPerSide() const=0;
    virtual int getNSidesPerCell() const=0;
    virtual int getNEdgesPerCell() const=0;

    virtual int getNNodes() const=0;
    virtual int getNCells() const=0;
    virtual int getNSides() const=0;
    virtual int getNEdges() const=0;

    virtual alat::Node getNodeOfCell(int iK) const=0;
    virtual alat::Node getNodeOfSide(int iS) const=0;

    virtual const Side getSide(int i) const=0;
    virtual const Cell getCell(int i) const=0;

    virtual void loadH5(const arma::hdf5_name& spec)=0;
    virtual void saveH5(const arma::hdf5_name& spec) const=0;
    virtual void readGmsh(std::string filename);
    virtual void writeVtk(std::string filename) const;
    virtual void writeBoundaryVtk(std::string filename) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

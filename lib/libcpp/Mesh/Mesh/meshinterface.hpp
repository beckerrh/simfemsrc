#ifndef __Mesh_MeshInterface_h
#define __Mesh_MeshInterface_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/map.hpp"
#include  "Mesh/enums.hpp"
#include  "Mesh/geometryobject.hpp"
#include  <memory>

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class BoundaryMeshUnit;
  class InterfaceMeshUnitsMap;
  class MeshUnitInterface;
}
namespace mesh
{
  class MeshInterface : public virtual alat::InterfaceBase
  {
  public:
    typedef alat::Map<int, std::shared_ptr<mesh::BoundaryMeshUnit> >  BoundaryMeshUnitsMap;

  protected:
    std::string getInterfaceName() const;

  public:
    ~MeshInterface();
    MeshInterface();
    MeshInterface( const MeshInterface& meshinterface);
    MeshInterface& operator=( const MeshInterface& meshinterface);

    virtual const MeshUnitInterface* getPlainMesh() const=0;
    virtual const MeshUnitInterface* getBoundaryMesh(int color) const=0;
    virtual const MeshUnitInterface* getInterfaceMesh(int color) const=0;
    virtual const BoundaryMeshUnitsMap& getBoundaryMeshUnitsMap() const=0;
    virtual const InterfaceMeshUnitsMap& getInterfaceMeshUnitsMap() const=0;
    virtual MeshUnitInterface* getPlainMesh()=0;
    virtual BoundaryMeshUnitsMap& getBoundaryMeshUnitsMap()=0;
    virtual InterfaceMeshUnitsMap& getInterfaceMeshUnitsMap()=0;

    virtual int getPartionId() const=0;
    virtual void setPartionId(int id)=0;
    virtual void addGeometryObject(meshEnums::geomobjtype type, const GeometryConstructorInterface* geometryconstructor=&trivialgeometryconstructor)=0;
    virtual std::string getInfo() const=0;
    virtual void exchangeInterfaces(int nproc)=0;

    virtual void loadH5(std::string filename)=0;
    virtual void saveH5(std::string filename) const=0;
    virtual void readGmsh(std::string filename);
    virtual void writeVtk(std::string filename) const;
    virtual void writeBoundaryVtk(std::string filename) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

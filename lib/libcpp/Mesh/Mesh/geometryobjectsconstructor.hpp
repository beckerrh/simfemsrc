#ifndef __Mesh_GeometryObjectsConstructor_h
#define __Mesh_GeometryObjectsConstructor_h

#include  "Alat/interfacebase.hpp"
#include  "Mesh/enums.hpp"
#include  <memory>

/*--------------------------------------------------------------------------*/

namespace mesh
{
  class GeometryObject;

  class GeometryObjectsConstructor : public alat::InterfaceBase
  {
  public:
    ~GeometryObjectsConstructor();
    GeometryObjectsConstructor();
    GeometryObjectsConstructor( const GeometryObjectsConstructor& geometryobjectsconstructor);
    GeometryObjectsConstructor& operator=( const GeometryObjectsConstructor& geometryobjectsconstructor);
    std::string getClassName() const;

    std::unique_ptr<mesh::GeometryObject> newGeometryObject(meshEnums::geomobjtype type) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

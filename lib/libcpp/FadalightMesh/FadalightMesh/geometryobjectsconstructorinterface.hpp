#ifndef __FadalightMesh_GeometryObjectsConstructorInterface_h
#define __FadalightMesh_GeometryObjectsConstructorInterface_h

#include  "Alat/interfacebase.hpp"
#include  <map>

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class GeometryObject;
  class GeometryObjectsConstructorInterface : public alat::InterfaceBase
  {
private:
protected:
    std::string getInterfaceName() const
    {
      return "GeometryObjectsConstructorInterface";
    }

public:
    ~GeometryObjectsConstructorInterface();
    GeometryObjectsConstructorInterface();
    GeometryObjectsConstructorInterface( const GeometryObjectsConstructorInterface& geometryobjectsconstructorinterface);
    GeometryObjectsConstructorInterface& operator=( const GeometryObjectsConstructorInterface& geometryobjectsconstructorinterface);
    virtual FadalightMesh::GeometryObject* newGeometryObject( const std::string& name ) const=0;
virtual void constructGeometryObject(std::map<std::string, FadalightMesh::GeometryObject*>& geo, const std::string& name) const = 0;
  };
}

/*--------------------------------------------------------------------------*/

#endif

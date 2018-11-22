#ifndef __FadalightMesh_GeometryObjectsConstructor_h
#define __FadalightMesh_GeometryObjectsConstructor_h

#include  "FadalightMesh/geometryobjectsconstructorinterface.hpp"
#include  "Alat/interfacebase.hpp"
#include  <map>
#include  <string>

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class GeometryObjectsConstructor : public virtual FadalightMesh::GeometryObjectsConstructorInterface
  {
public:
	GeometryObjectsConstructor();
    std::string getClassName() const;
FadalightMesh::GeometryObject* newGeometryObject( const std::string& name ) const;
   void constructGeometryObject(std::map<std::string, FadalightMesh::GeometryObject*>& geo, const std::string& name) const;
   };
}

/*--------------------------------------------------------------------------*/

#endif

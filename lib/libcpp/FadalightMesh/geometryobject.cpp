#include  "FadalightMesh/geometryobject.hpp"
#include  <cassert>

using namespace FadalightMesh;

/*--------------------------------------------------------------------------*/
GeometryObject::~GeometryObject(){}
GeometryObject::GeometryObject() : alat::InterfaceBase(){}
GeometryObject::GeometryObject( const GeometryObject& geometryobject) : alat::InterfaceBase(geometryobject)
{
  assert(0);
}
GeometryObject& GeometryObject::operator=( const GeometryObject& geometryobject)
{
  InterfaceBase::operator=(geometryobject);
  assert(0);
  return *this;
}
std::string GeometryObject::getInterfaceName() const
{
  return "GeometryObject";
}

/*--------------------------------------------------------------------------*/

void GeometryObject::constructGeometryObject(const FadalightMesh::MeshInterface* mesh)
{
  _notWritten("constructGeometryObject");
}

#include  "Mesh/geometryobject.hpp"
#include  <cassert>

using namespace mesh;

/*--------------------------------------------------------------------------*/
GeometryObject::~GeometryObject(){}
GeometryObject::GeometryObject() : alat::InterfaceBase()
{
  // std::cerr << "GeometryObject::GeometryObject()\n";
}
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
void GeometryObject::construct(const mesh::MeshUnitInterface* mesh, const GeometryConstructorInterface* geometryconstructor)
{
  _notWritten("construct");
}

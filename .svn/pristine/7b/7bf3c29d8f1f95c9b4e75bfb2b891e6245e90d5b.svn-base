#include  "Alat/map.hpp"
#include  "Mesh/geometryobjectsconstructor.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Mesh/measureofcell.hpp"
#include  "Mesh/normals.hpp"
#include  "Mesh/nodescellsweight.hpp"
#include  <string>

using namespace mesh;

/*--------------------------------------------------------------------------*/
GeometryObjectsConstructor::~GeometryObjectsConstructor() {}
GeometryObjectsConstructor::GeometryObjectsConstructor() : alat::InterfaceBase() {}
GeometryObjectsConstructor::GeometryObjectsConstructor(const GeometryObjectsConstructor& geometryobjectsconstructor) : alat::InterfaceBase(geometryobjectsconstructor) {}
GeometryObjectsConstructor& GeometryObjectsConstructor::operator=( const GeometryObjectsConstructor& geometryobjectsconstructor)
{
  assert(0);
  return *this;
}
std::string GeometryObjectsConstructor::getClassName() const
{
  return "GeometryObjectsConstructor";
}

/*--------------------------------------------------------------------------*/
std::unique_ptr<mesh::GeometryObject> GeometryObjectsConstructor::newGeometryObject(meshEnums::geomobjtype type) const
{
  // std::cerr << "GeometryObjectsConstructor::newGeometryObject()=\n";
  if(type==meshEnums::MeasureOfCell)
  {
    return std::unique_ptr<mesh::GeometryObject>(new MeasureOfCell);
  }
  else if(type==meshEnums::Normals)
  {
    return std::unique_ptr<mesh::GeometryObject>(new Normals);
  }
  else if(type==meshEnums::NodesCellsWeight)
  {
    return std::unique_ptr<mesh::GeometryObject>(new NodesCellsWeight);
  }
  else if(type==meshEnums::CutInterface)
  {
    return std::unique_ptr<mesh::GeometryObject>(new CutInterface);
  }
  _error_string("newGeometryObject", "unknown GeometryObject", meshEnums::geomObjTypeToString(type));
  return std::unique_ptr<mesh::GeometryObject>(nullptr);
}

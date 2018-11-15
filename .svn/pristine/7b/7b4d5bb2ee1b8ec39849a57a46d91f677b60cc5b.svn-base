#include  "FadalightMesh/coarseninfo.hpp"
#include  "FadalightMesh/curvedboundaryinformation.hpp"
#include  "FadalightMesh/geometryobjectsconstructor.hpp"
#include  "FadalightMesh/hangingnodeinfo.hpp"
#include  "FadalightMesh/hangingsideinfo.hpp"
#include  "FadalightMesh/refineinfo.hpp"
#include  <string>

using namespace FadalightMesh;

 /*--------------------------------------------------------------------------*/
 GeometryObjectsConstructor::GeometryObjectsConstructor() : FadalightMesh::GeometryObjectsConstructorInterface() {}
 std::string GeometryObjectsConstructor::getClassName() const
    {
      return "GeometryObjectsConstructor";
    }

/*--------------------------------------------------------------------------*/
void FadalightMesh::GeometryObjectsConstructor::constructGeometryObject(std::map<std::string, FadalightMesh::GeometryObject*>& geo, const std::string& name ) const
{
  geo[name] = newGeometryObject(name);
}

/*--------------------------------------------------------------------------*/
FadalightMesh::GeometryObject* GeometryObjectsConstructor::newGeometryObject( const std::string& name ) const
{

  if(name == "HangingSideInfo")
  {
    return  new FadalightMesh::HangingSideInfo;
  }
  else if(name == "HangingNodeInfo")
  {
    return new FadalightMesh::HangingNodeInfo;
  }
  else if(name == "RefineInfo")
  {
    return  new FadalightMesh::RefineInfo;
  }
  else if(name == "CoarsenInfo")
  {
    return  new FadalightMesh::CoarsenInfo;
  }
  else if(name == "CurvedBoundaryInformation")
  {
    return new FadalightMesh::CurvedBoundaryInformation;
  }
    _notWritten("createGeometryObject("+name+")");
}

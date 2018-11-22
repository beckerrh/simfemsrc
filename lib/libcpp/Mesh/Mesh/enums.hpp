#ifndef  __Mesh_enum_h
#define  __Mesh_enum_h

#include  <string>

/*---------------------------------------------*/
namespace meshEnums
{
  enum meshtype {PointMesh, LineMesh, TriangleMesh, QuadrilateralMesh, HexahedralMesh, TetrahedralMesh};
  std::string meshTypeToString(meshtype m);
  meshtype stringTomeshType(std::string s);

  enum geomobjtype {MeasureOfCell, Normals, NodesCellsWeight, CutInterface};
  std::string geomObjTypeToString(geomobjtype m);
  geomobjtype stringTogeomObjType(std::string s);

  enum meshunittype {None, Plain, Boundary, Interface};
  std::string meshunitTypeToString(meshunittype m);
}

#endif

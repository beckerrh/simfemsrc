#ifndef  __FadalightMesh_enum_h
#define  __FadalightMesh_enum_h

#include  <string>

/*---------------------------------------------*/

namespace FadalightMeshEnums
{
  enum meshtype {LineMesh, TriangleMesh, QuadrilateralMesh, HexahedralMesh, TetrahedralMesh, CouplingMesh2D,CouplingMesh3D};

  std::string meshTypeToString(meshtype m);
}

#endif

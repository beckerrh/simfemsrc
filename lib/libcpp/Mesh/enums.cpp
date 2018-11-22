#include  "Mesh/enums.hpp"
#include  <assert.h>
#include  <iostream>

namespace meshEnums
{
/*---------------------------------------------*/
  std::string meshTypeToString(meshtype m)
  {
    if(m == PointMesh) {return "PointMesh";}
    else if(m == LineMesh) {return "LineMesh";}
    else if(m == TriangleMesh) {return "TriangleMesh";}
    else if(m == QuadrilateralMesh) {return "QuadrilateralMesh";}
    else if(m == HexahedralMesh) {return "HexahedralMesh";}
    else if(m == TetrahedralMesh) {return "TetrahedralMesh";}
    return "none";
  }
  meshtype stringTomeshType(std::string s)
  {
    if(s == "PointMesh") {return PointMesh;}
    else if(s == "LineMesh") {return LineMesh;}
    else if(s == "TriangleMesh") {return TriangleMesh;}
    else if(s == "QuadrilateralMesh") {return QuadrilateralMesh;}
    else if(s == "HexahedralMesh") {return HexahedralMesh;}
    else if(s == "TetrahedralMesh") {return TetrahedralMesh;}
    std::cerr<<"****stringTomeshType: meshtype not defined "<<s<<'\n';
    assert(0);
    return TetrahedralMesh;
  }

/*---------------------------------------------*/
  std::string geomObjTypeToString(geomobjtype m)
  {
    if(m == MeasureOfCell) {return "MeasureOfCell";}
    else if(m == Normals) {return "Normals";}
    else if(m == NodesCellsWeight) {return "NodesCellsWeight";}
    else if(m == CutInterface) {return "CutInterface";}
    return "none";
  }
  geomobjtype stringTogeomObjType(std::string s)
  {
    if(s == "MeasureOfCell") {return MeasureOfCell;}
    else if(s == "Normals") {return Normals;}
    else if(s == "NodesCellsWeight") {return NodesCellsWeight;}
    else if(s == "CutInterface") {return CutInterface;}
    std::cerr<<"****stringTogeomObjType: geomobjtype not defined "<<s<<'\n';
    assert(0);
    return MeasureOfCell;
  }
  /*---------------------------------------------*/
  std::string meshunitTypeToString(meshunittype m)
  {
    if(m == None) {return "None";}
    else if(m == Plain) {return "Plain";}
    else if(m == Boundary) {return "Boundary";}
    else if(m == Interface) {return "Interface";}
    return "none";
  }
}

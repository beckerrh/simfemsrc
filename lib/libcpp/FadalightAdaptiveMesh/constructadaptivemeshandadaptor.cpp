#include  "FadalightAdaptiveMesh/basicadaptivemesh2d.hpp"
// #include  "FadalightAdaptiveMesh/basicadaptivemesh3d.hpp"
#include  "FadalightAdaptiveMesh/coarsener2d.hpp"
// #include  "FadalightAdaptiveMesh/coarsener3d.hpp"
#include  "FadalightAdaptiveMesh/constructadaptivemeshandadaptor.hpp"
#include  "FadalightAdaptiveMesh/refiner2d.hpp"
// #include  "FadalightAdaptiveMesh/refiner3d.hpp"
#include  "FadalightMesh/getmeshtype.hpp"

using namespace FadalightAdaptiveMesh;
using namespace std;

/*--------------------------------------------------------------------------*/

ConstructAdaptiveMeshAndAdaptor::ConstructAdaptiveMeshAndAdaptor(AdaptiveMeshInterface*& M, MeshAdaptorInterface*& AM, const std::string& meshname, const std::string& adaption_type)
{
  pair<string,string> p = FadalightMesh::getMeshType(meshname);
  string type = p.first;
  string datatype = p.second;

  if(type == "FadalightMesh::TriangleMesh")
  {
    M = new BasicAdaptiveMesh2d<3>;
  }
  else if(type == "FadalightMesh::QuadrilateralMesh")
  {
    M = new BasicAdaptiveMesh2d<4>;
  }
  else if(type == "FadalightMesh::HexahedralMesh")
  {
    assert(0);
    // M = new BasicAdaptiveMesh3d<8,6,12,4>;
  }
  else
  {
    std::cerr<<"***Error: ConstructAdaptiveMeshAndAdaptor: bad  type: "<<type<<'\n';
    assert(0);
    exit(1);
  }
  if(type == "FadalightMesh::HexahedralMesh")
  {
    assert(0);
    // if (adaption_type == "refine")
    // {
    //    AM = new Refiner3d(M);
    // }
    // else if(adaption_type == "coarse")
    // {
    //   AM = new  Coarsener3d(M);
    // }
    // else
    // {
    //   std::cerr<<"***Error: ConstructAdaptiveMeshAndAdaptor: bad adaption type: "<<type<<"  "<<adaption_type<<'\n';
    //   assert(0);
    //   exit(1);
    // }
  }
  else if ((type == "FadalightMesh::QuadrilateralMesh") or (type == "FadalightMesh::TriangleMesh"))
  {
     if(adaption_type == "coarse")
     {
       AM = new  Coarsener2d(M);
     }
     else if(adaption_type == "refine")
     {
       AM = new  Refiner2d(M);
     }
     else
     {
       std::cerr<<"***Error: ConstructAdaptiveMeshAndAdaptor: bad adaption type: "<<adaption_type<<'\n';
       assert(0);
       exit(1);
     }
  }
}

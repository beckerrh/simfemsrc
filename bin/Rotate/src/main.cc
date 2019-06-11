#include  "Alat/map.hpp"
#include  "Alat/set.hpp"
#include  "Alat/tokenize.hpp"
#include  "FadalightMesh/hexahedralmesh.hpp"
#include  "FadalightMesh/quadrilateralmesh.hpp"
#include  <algorithm>
#include  <iostream>
#include  <math.h>

/*--------------------------------------------------------------------------*/

int main(int argc, char** argv)
{
  if(argc != 6)
  {
    std::cerr << "*** ERROR in " << argv[0] <<" usage: filein  fileout meshtype datatype angle\n";
    return 1;
  }
  std::string infilename(argv[1]);
  std::string outfilename(argv[2]);
  std::string mesh_type(argv[3]);
  std::string datatype(argv[4]);
  double angle = atof(argv[5]);

  FadalightMesh::MeshInterface* mesh;
  if(mesh_type == "quad")
  {
    mesh = new FadalightMesh::QuadrilateralMesh;
    FadalightMesh::QuadrilateralMesh* quadmesh = dynamic_cast<FadalightMesh::QuadrilateralMesh*>( mesh );
    infilename += ".quad";
    quadmesh->readQuad(infilename);
  }
  else if(mesh_type == "hex")
  {
    mesh = new FadalightMesh::HexahedralMesh;
    FadalightMesh::HexahedralMesh* hexmesh = dynamic_cast<FadalightMesh::HexahedralMesh*>( mesh );
    infilename += ".hex";
    hexmesh->readHex(infilename);
  }
  else
  {
    std::cerr<<"***ERROR in " << argv[0] <<"  : invalid mesh type. Must be quad, hex or tri, but "<<mesh_type<<" given."<<'\n';
    assert(0);
    exit(1);
  }

  alat::Vector<alat::Node>& nodes = mesh->getAllNodes();
  
  double c = cos(2.0*M_PI*angle/360.0);
  double s = sin(2.0*M_PI*angle/360.0);
  for(int i=0;i<nodes.size();i++)
  {
    alat::Node& node = nodes[i];
    double x = node.x();
    double y = node.y();
    double z = node.z();
    node.x() = c*x-s*y;
    node.y() = s*x+c*y;
  }


  mesh->writeVtk(outfilename);
  mesh->writeBoundaryVtk(outfilename);
  mesh->write(outfilename);

  delete mesh;

  return 0;
}

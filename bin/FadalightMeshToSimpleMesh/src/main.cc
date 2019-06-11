#include  "FadalightMesh/hexahedralmesh.hpp"
#include  "FadalightMesh/quadrilateralmesh.hpp"
// #include  "FadalightMesh/trianglemesh.hpp"
#include  <iostream>

/* ---------------------------------- */

int main(int argc, char** argv)
{
  if(argc != 5)
  {
    std::cerr << "*** ERROR in " << argv[0] <<" usage: filein  fileout meshtype datatype\n";
    return 1;
  }
  std::string infilename(argv[1]);
  std::string outfilename(argv[2]);
  std::string mesh_type(argv[3]);
  std::string datatype(argv[4]);
  FadalightMesh::MeshInterface* mesh;
  if(mesh_type == "quad")
  {
    mesh = new FadalightMesh::QuadrilateralMesh;
    outfilename += ".quad";
  }
  else if(mesh_type == "hex")
  {
    mesh = new FadalightMesh::HexahedralMesh;
    outfilename += ".hex";
  }
  else
  {
    std::cerr<<"***ERROR in " << argv[0] <<"  : invalid mesh type. Must be quad, hex or tri, but "<<mesh_type<<" given."<<'\n';
    assert(0);
    exit(1);
  }
  mesh->readFadalightMesh(infilename);
  mesh->write(outfilename);
  // std::cerr << "@@@@@@@@@@@@@@@@@@@@ " << argv[0] <<  ": " << infilename << " ====>>> " << outfilename << "\n";
  // assert(0);
  delete mesh;
}

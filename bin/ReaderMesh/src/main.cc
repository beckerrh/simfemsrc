#include  "FadalightMesh/hexahedralmesh.hpp"
#include  "FadalightMesh/quadrilateralmesh.hpp"
#include  "FadalightMesh/quadtotri.hpp"
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
  FadalightMesh::MeshInterface* M;
  if(mesh_type == "quad")
  {
    M = new FadalightMesh::QuadrilateralMesh;
    FadalightMesh::QuadrilateralMesh* Mquad = dynamic_cast<FadalightMesh::QuadrilateralMesh*>( M );
    infilename += ".quad";
    Mquad->readQuad(infilename);
  }
  else if(mesh_type == "hex")
  {
    M = new FadalightMesh::HexahedralMesh;
    FadalightMesh::HexahedralMesh* Mhex = dynamic_cast<FadalightMesh::HexahedralMesh*>( M );
    infilename += ".hex";
    Mhex->readHex(infilename);
  }
  else if(mesh_type == "quadtotri")
  {
    M = new FadalightMesh::QuadToTri(infilename);
  }
  // else if(mesh_type == "tri")
  // {
  //   M = new FadalightMesh::TriangleMesh;
  //   FadalightMesh::TriangleMesh* Mtri = dynamic_cast<FadalightMesh::TriangleMesh*>( M );
  //   Mtri->ReadTri(infilename);
  // }
  else
  {
    std::cerr<<"***ERROR in " << argv[0] <<"  : invalid mesh type. Must be quad, hex or tri, but "<<mesh_type<<" given."<<'\n';
    assert(0);
    exit(1);
  }
  if(datatype == "ascii")
  {
    M->writeFadalightMesh(outfilename, arma::arma_ascii);
    M->writeVtk(outfilename);
  }
  else
  {
    M->writeFadalightMesh(outfilename);    
  }
  M->writeBoundaryVtk(outfilename);

  delete M;

  return 0;
}

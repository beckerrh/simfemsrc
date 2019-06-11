#include  "FadalightMesh/refineinfo.hpp"
#include  "FadalightMesh/quadtotri.hpp"
#include  "FadalightMesh/multilevelmesh.hpp"
#include  "FadalightMesh/getmeshtype.hpp"
#include  <sstream>

/*----------------------------------------------------------*/

int main(int argc, char** argv)
{

  if( argc != 3)
  {
    std::cerr << "*** " << argv[0] << " usage : trimesh (FadalightMesh) datatype\n";
    return 1;
  }
  std::string trimesh(argv[1]);
  std::string datatype(argv[2]);

  FadalightMesh::QuadToTri quadtotri;
  quadtotri.readFadalightMesh(trimesh);
  FadalightMesh::MultiLevelMesh mlmeshquad("FadalightMesh::QuadrilateralMesh");
  mlmeshquad.readFadalightMesh(quadtotri.getFileNameQuad());
  int nlevels = mlmeshquad.getNLevels();
  
  std::string multimeshname = trimesh+".fadalightmesh/MultiMesh";
  std::string cmd = "rm -rf "+multimeshname;
  int error = system( cmd.c_str() );
  assert(!error);
  cmd = "mkdir "+multimeshname;
  error = system( cmd.c_str() );
  assert(!error);
  std::string filename = multimeshname+"/n";
  std::ofstream file( filename.c_str() );
  assert( file.is_open() );
  file<<nlevels<<"\n";
  file.close();
  
  
  for(int level=1;level<nlevels;level++)
  {
    std::string quadmeshname = mlmeshquad.getMesh(level)->getInfileName();
    // std::cerr << "quadmeshname = "<< quadmeshname << "\n";
    FadalightMesh::QuadToTri trimesh =  FadalightMesh::QuadToTri(quadmeshname);
    alat::StringPair p = FadalightMesh::getMeshType(quadmeshname);
    
    // std::cerr << "??????????????????????????????? reste a construire le refinfo !!!!!!!!!!\n";
    std::stringstream ss0;
    ss0<<level;
    std::string outfilename = multimeshname+"/mesh."+ss0.str();
    if(datatype == "ascii")
    {
      trimesh.writeVtk(outfilename);
      trimesh.writeFadalightMesh(outfilename, arma::arma_ascii);
    }
    else
    {
      trimesh.writeFadalightMesh(outfilename);      
    }
  }
  return 0;
}

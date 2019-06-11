#include  "FadalightMesh/quadrilateralmesh.hpp"
#include  "FadalightMesh/trianglemesh.hpp"
#include  <fstream>
#include  "Alat/set.hpp"

using namespace std;

/*----------------------------------------------------------*/

int main(int argc, char** argv)
{   
  if(argc != 3)
  {
    std::cerr << "*** " << argv[0] << " usage : infilename outfilename (FadalightMesh)\n";
    return 1;
  }   
  std::string infilename(argv[1]);
  std::string outfilename(argv[2]);

  FadalightMesh::QuadrilateralMesh M;
  M.readFadalightMesh(outfilename);

  // read cells quad --> tri
  alat::armaivec quad_of_tri;

  std::string filename = infilename + ".fadalightmesh/quad_of_tri";
  std::cerr<<"filename= "<<filename<<std::endl;
  std::ifstream file(filename.c_str()); 
  assert(file.is_open());
  quad_of_tri.load(file);
  file.close();   
  
  // read marked cells
  alat::armaivec marked_cells;
  filename = infilename + ".fadalightmesh/marked_cells";
  file.open(filename.c_str());
  if(not file.is_open())
  {
    std::cerr << "*** ERROR in " << argv[0] << ": cannot open file \"" << filename << "\"\n";
    assert(0);
    exit(1);
  }
  marked_cells.load(file);
  file.close();
  
  alat::IntSet new_marked;
  for(int i=0;i<marked_cells.size();i++)
  {    
    new_marked.insert(quad_of_tri[marked_cells[i]]);
  }
  std::cout << "*** marked (quad) : " << new_marked.size() << "\n";
  
  alat::armaivec marked_cells_quad(new_marked.size());                                                       
  int count=0;
  for(alat::IntSet::const_iterator p=new_marked.begin();p!=new_marked.end();p++) marked_cells_quad[count++]=*p;

  // write marked cells
  {
    filename = outfilename + ".fadalightmesh/marked_cells";
    std::ofstream file(filename.c_str());
    assert(file.is_open());
    // marked_cells_quad.save(file,"ascii");
    marked_cells_quad.save(file,arma::arma_ascii);
    file.close();
  }

  return 0;
}

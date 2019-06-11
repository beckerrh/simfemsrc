// #include  "Alat/iomanager.hpp"
// #include  "Fadalight/solutionio.hpp"
#include  "Alat/tokenize.hpp"
// #include  "Liniowa/variablevector.hpp"
#include  "FadalightMesh/markercoarsepercentage.hpp"
#include  "FadalightMesh/markermean.hpp"
#include  "FadalightMesh/markeroptimal.hpp"
#include  "FadalightMesh/markeroptimalcoarse.hpp"
#include  "FadalightMesh/markerpercentage.hpp"
#include  "FadalightMesh/markerrandom.hpp"
#include  <dirent.h>
#include  <stdlib.h>
#include  <fstream>

using namespace std;

/// This program generates a set of marked cells from cell-wise values of an error estimator.
/// For this purpose we suppose that the run-directory contains a file 'Estimator.fadalightsolution'.
/// It produces as output a file 'marked_cells' in the mesh-directory

/*----------------------------------------------------------*/

int getNumberOfCells(std::string inmeshname)
{
  std::string meshname = inmeshname+".fadalightmesh/name";
  ifstream meshnamefile( meshname.c_str() );
  if( not meshnamefile.is_open() )
  {
    std::cerr<<"*** ERROR in getNumberOfCells(): cannot open file \""<<meshname<<"\"\n";
    assert(0);
    exit(1);
  }
  std::string toto;
  meshnamefile>>toto;
  meshnamefile>>toto;
  meshnamefile>>toto;
  vector<string> bouts = alat::Tokenize(toto, "_");
  assert(bouts.size() > 1);
  meshnamefile.close();
  return atoi( bouts[bouts.size()-2].c_str() );
}

/*----------------------------------------------------------*/

int getDimension(std::string inmeshname)
{
  std::string meshname = inmeshname+".fadalightmesh/name";
  ifstream meshnamefile( meshname.c_str() );
  if( not meshnamefile.is_open() )
  {
    std::cerr<<"*** ERROR in getNumberOfCells(): cannot open file \""<<meshname<<"\"\n";
    assert(0);
    exit(1);
  }
  std::string toto;
  meshnamefile>>toto;
  meshnamefile.close();
  int dimension = 3;
  if( ( toto == "FadalightMesh::QuadrilateralMesh" )||( toto == "FadalightMesh::TriangleMesh" ) )
  {
    dimension = 2;
  }
  else if(toto == "FadalightMesh::HexaedralMesh" )
  {
    dimension = 3;
  }
  else
  {
    std::cerr << "*** ERROR in getDimension(): unknown mesh " << toto << "\n";
    assert(0);
    exit(1);
  }
  return dimension;
}

// /*----------------------------------------------------------*/
//
// void readFromFile(Liniowa::VariableVector& estimator, std::string filename)
// {
//     assert(0);
//     exit(1);
// //  Fadalight::SolutionDescription _soldesc;
// //  _soldesc.read(filename+".desc");
// //  int ncomp = _soldesc.getNComponents();
// //  estimator.setNComponents( ncomp );
// //  estimator.set_size( _soldesc.getNDofPerComponent() );
// //  std::ifstream file( filename.c_str() );
// //  std::cerr << "readFromFile() faudrait faire nouvelle version\n";
// //  assert(0);
// //  // estimator.readBlock( file, 0, ncomp, _soldesc.getDataType(), _soldesc.getPrecision() );
// //  file.close();
// }

/*----------------------------------------------------------*/

int main(int argc, char** argv)
{
  if(argc < 6)
  {
    std::cerr << "***ERROR in "<<argv[0]<<" expected parameters : inmeshname estimatorfile datatype marker marking_parameter\n";
    std::cerr << "I got " << argc << "\n";
    std::cerr << "I got " << argv << "\n";
    return 1;
  }
  std::string inmeshname(argv[1]);
  std::string estimatorfile(argv[2]);
  std::string datatype(argv[3]);
  std::string marker(argv[4]);
  double marking_parameter( atof(argv[5]) );
  // std::cerr << "\n\n XXXXXXXXXXXXXXXXXXX marker " << marking_parameter << "\n";
  

  FadalightMesh::MarkerInterface* MARKER;
  // std::cerr << "\n\n XXXXXXXXXXXXXXXXXXX marker " << marker << "\n";
  if(marker == "percentage")
  {
    MARKER = new FadalightMesh::MarkerPercentage(marking_parameter);
  }
  else if(marker == "coarsepercentage")
  {
    MARKER = new FadalightMesh::MarkerCoarsePercentage(marking_parameter);
  }
  else if(marker == "mean")
  {
    MARKER = new FadalightMesh::MarkerMean(marking_parameter);
  }
  else if(marker == "optimal")
  {
    int dimension = getDimension(inmeshname);
    MARKER = new FadalightMesh::MarkerOptimal(marking_parameter, dimension);
  }
  else if(marker == "optimalcoarse")
  {
    int dimension = getDimension(inmeshname);
    MARKER = new FadalightMesh::MarkerOptimalCoarse(marking_parameter, dimension);
  }
  else if(marker == "random")
  {
    int n = getNumberOfCells(inmeshname);
    MARKER = new FadalightMesh::MarkerRandom(n, marking_parameter);
  }
  else
  {
    std::cerr << "*** ERROR in " <<argv[0] << " unknown marker \"" << marker << "\"\n";
    assert(0);
    exit(1);
  }
  std::string outfilename = inmeshname+".fadalightmesh/marked_cells";
  if(datatype == "ascii")
  {
    MARKER->write(outfilename, arma::arma_ascii);    
  }
  else
  {
    MARKER->write(outfilename);    
  }
  delete MARKER;
}

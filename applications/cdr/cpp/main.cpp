#include  <Mesh/mesh.hpp>
#include  <Mesh/enums.hpp>
#include  "solver.hpp"
#include  <Solvers/options.hpp>
#include  <list>

/*---------------------------------------------------------*/
void test(std::string testname, meshEnums::meshtype type, std::string method, std::string fem)
{
  double diff=1.0, alpha=0.0, deltasupg=0.5;
  // std::string beta = "zero";
  std::string beta = "one";
  if(testname=="poisson")
  {
  }
  else if(testname=="cdrexp")
  {
    std::string beta = "east";
    diff = 0.001;
  }
  else if(testname=="rdcosh")
  {
    alpha = 1.0;
    diff = 0.0001;
    deltasupg = 0.5;
  }
  else
  {
    std::cerr << "unknwown testname "<<testname << "\n";
  }
  // std::string application = "constant";
  std::string application = "linear";
  // std::string application = "quadratic";
  // std::string application = "cdexplayer";
  // std::string application = "rdcosh";

  int partion_id=1;
  bool construct_pdepartmeshes=true;
  std::shared_ptr<mesh::MeshInterface> mesh = mesh::Mesh::create(type, partion_id, construct_pdepartmeshes);
  std::string meshfile;
  if(type==meshEnums::LineMesh)
  {
    meshfile = "unitline";
  }
  else if(type==meshEnums::TriangleMesh)
  {
    meshfile = "unitsquare";
  }
  else if(type==meshEnums::TetrahedralMesh)
  {
    meshfile = "unitcube";
  }
  else
  {
    std::cerr << "unknown meshtype\n";
    exit(1);
  }
  mesh->readGmsh("/Users/becker/Programs/simfem/installdir/meshes/"+meshfile+".msh");
  mesh->addGeometryObject(meshEnums::MeasureOfCell);
  mesh->addGeometryObject(meshEnums::Normals);
  mesh->addGeometryObject(meshEnums::NodesCellsWeight);
  std::cerr << "Mesh:\n" << mesh->getInfo() << "\n";

  Solver solver;
  solver.setMesh(mesh);
  solver.setParameter("method", method);
  solver.setParameter("application", application);
  solver.setParameter("beta", beta);
  solver.setParameter("fem", fem);
  solver.setParameter("diff", diff);
  solver.setParameter("alpha", alpha);
  solver.setParameter("deltasupg", deltasupg);
  solver.setOutputOptions(solver_options::datadir) = meshfile;
  solver.init();
  std::cerr << "Solver:\n" << solver.getInfo() << "\n";
  solvers::ErrorsMap errors =solver.run();
  std::cerr << "Errors:\n" << errors << "\n";
  solver.writeXdmf();
}

/*---------------------------------------------------------*/
int main(int argc, char** argv)
{
  std::string fem = "P1";
  // std::string fem = "CR1";
  // std::string fem = "P2";
  std::string testname="poisson";
  std::list<std::string> pdeparts;
  // pdeparts.push_back("traditional");
  // pdeparts.push_back("traditionalintegration");

  // pdeparts.push_back("nitsche");
  // pdeparts.push_back("nitscheintegration");
  pdeparts.push_back("newnitscheintegration");
  // pdeparts.push_back("traditionalintegration");
  for(std::list<std::string>::const_iterator p=pdeparts.begin();p!=pdeparts.end();p++)
  {
    std::cerr << "##############  " << *p <<  "  ##############\n";
    test(testname, meshEnums::LineMesh, *p, fem);
    test(testname, meshEnums::TriangleMesh, *p, fem);
    test(testname, meshEnums::TetrahedralMesh, *p, fem);
  }
  return 0;
}

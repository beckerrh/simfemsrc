#include  <Mesh/mesh.hpp>
#include  <Mesh/enums.hpp>
#include  <Solvers/solver.hpp>
#include  <Solvers/options.hpp>

/*---------------------------------------------------------*/
void test(meshEnums::meshtype type, bool armamatrix)
{
  int partion_id=1;
  bool construct_bdrymeshes=true;
  std::shared_ptr<mesh::MeshInterface> mesh = mesh::Mesh::create(type, partion_id, construct_bdrymeshes);
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
  mesh->readGmsh("@simfeminstallpath@/meshes/"+meshfile+".msh");
  mesh->addGeometryObject(meshEnums::MeasureOfCell);
  mesh->addGeometryObject(meshEnums::Normals);
  mesh->addGeometryObject(meshEnums::NodesCellsWeight);
  // mesh->saveH5("toto.h5");
  std::cerr << "Mesh:\n" << mesh->getInfo() << "\n";
  // mesh->loadH5("mesh.h5");

  solver_options::opts opts(solver_options::flag_none);// = solver_options::strongdir;
  if(armamatrix) opts += solver_options::armamat;
  solvers::Solver solver(mesh, opts);
  solver.setOutputOptions(solver_options::datadir) = meshfile;
  // solver.addVariablePlain("U", 1, solverEnums::fem::P1);
  solver.init();
  std::cerr << "Solver:\n" << solver.getInfo() << "\n";
  solver.run();
  solver.saveSolution();
  solver.writeXdmf();
}

/*---------------------------------------------------------*/
int main(int argc, char** argv)
{
  bool armamatrix = false;
  if(argc==2)
  {
    armamatrix = atoi(argv[1]);
  }

  test(meshEnums::LineMesh, armamatrix);
  test(meshEnums::TriangleMesh, armamatrix);
  test(meshEnums::TetrahedralMesh, armamatrix);
  return 0;
}

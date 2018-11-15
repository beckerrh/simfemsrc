#include  <Mesh/mesh.hpp>
#include  <Mesh/enums.hpp>
#include  "solver.hpp"
#include  <Solvers/options.hpp>

/*---------------------------------------------------------*/
void test(meshEnums::meshtype type)
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
  std::cerr << "Mesh:\n" << mesh->getInfo() << "\n";

  solver_options::opts opts = solver_options::dynamic;
  Solver solver(mesh, opts);
  solver.setOutputOptions(solver_options::datadir) = meshfile;
  mesh::TimeMeshData timemeshdata;
  timemeshdata.n = 1;
  timemeshdata.last = 40.;
  solver.setTimeMesh(timemeshdata);
  solver.init();
  std::cerr << "Solver: " << solver.getInfo() << "\n";
  solver.run();
  solver.writeXdmf();
}

/*---------------------------------------------------------*/
int main(int argc, char** argv)
{
  // test(meshEnums::LineMesh);
  test(meshEnums::TriangleMesh);
  // test(meshEnums::TetrahedralMesh);
  return 0;
}

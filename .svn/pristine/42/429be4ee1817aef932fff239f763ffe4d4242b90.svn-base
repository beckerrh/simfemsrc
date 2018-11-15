#include  <Mesh/mesh.hpp>
#include  <Mesh/enums.hpp>
#include  "solver.hpp"
#include  <Solvers/options.hpp>

/*---------------------------------------------------------*/
void test(std::string method)
{
  meshEnums::meshtype type=meshEnums::TriangleMesh;
  // meshEnums::meshtype type=meshEnums::TetrahedralMesh;
  int partion_id=1;
  std::shared_ptr<mesh::MeshInterface> mesh = mesh::Mesh::create(type, partion_id);
  std::string meshfile;
  if(type==meshEnums::TriangleMesh) {meshfile = "unitsquare";}
  else if(type==meshEnums::TetrahedralMesh) {meshfile = "unitcube";}
  mesh->readGmsh("@simfeminstallpath@/meshes/"+meshfile+".msh");
  mesh->addGeometryObject(meshEnums::MeasureOfCell);
  mesh->addGeometryObject(meshEnums::Normals);

  Solver solver(mesh);
  solver.setOutputOptions(solver_options::datadir) = meshfile;

  // solver.setParameter("fem", "cr1");
  solver.setParameter("kin", 1.0);
  // solver.setParameter("kex", 2.0);
  solver.setParameter("kex", 1e6);
  // solver.setParameter("kex", 1.0);
  solver.setParameter("xgamma", 0.5);
  solver.setParameter("gamma", 1.0);
  solver.setParameter("application", "linear_straight");
  // solver.setParameter("application", "linear");
  // solver.setParameter("application", "quadratic_circle");
  solver.setParameter("application", "cardioide");
  solver.setParameter("method", method);

  solver.init();
  std::cerr<< solver.getInfo();
  solvers::ErrorsMap errors =solver.run();
  std::cerr << "Errors:\n" << errors << "\n";
  solver.writeXdmf();
}

/*---------------------------------------------------------*/
int main(int argc, char** argv)
{
  // test("nitsche");
  // test("nitschestabc");
  // test("nitschestabnc");
  // test("newnitsche");
  // test("newnitsche2");
  // test("strong");
  test("P1IFEM");
  // test("CR1IFEM");
  return 0;
}

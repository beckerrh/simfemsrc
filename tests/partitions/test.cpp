#include  <mpi.h>
#include  <iomanip>
#include  <sstream>
#include  <Mesh/mesh.hpp>
#include  <Solvers/solver.hpp>

/*---------------------------------------------------------*/
int main(int argc, char** argv)
{
  int nproc, iproc, name_len;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  bool mpirun = true;
  if(argc>1 and std::string(argv[1])=="nompi")
  {
    mpirun=false;
    nproc = 2;
    iproc = 1;
    if(argc==4)
    {
      iproc = atoi(argv[2]);
      nproc = atoi(argv[3]);
    }
  }

  if(mpirun)
  {
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
    MPI_Get_processor_name(processor_name, &name_len);
  }

  if(iproc==0)
  {
    std::string mesh = "unitsquare";
    std::string geofile = "@simfeminstallpath@/meshes/"+mesh+".geo";
    std::stringstream cmd;
    cmd << "gmsh -2  " << geofile << " -o " << mesh<< ".msh -part " <<nproc-1 << " -oneFilePerPart";
    system(cmd.str().c_str());
    std::cerr << "mesh ready\n";
  }
  if(mpirun) {MPI_Barrier(MPI_COMM_WORLD);}

  std::cerr << "I am " << iproc << " of " << nproc << " name is " << processor_name << "\n";

  std::shared_ptr<mesh::MeshInterface> mesh;
  if(iproc)
  {
    std::stringstream ss;
    ss << "unitsquare.msh_";
    ss << std::setfill('0') << std::setw(6) << iproc;
    meshEnums::meshtype type = meshEnums::TriangleMesh;
    bool construct_bdrymeshes = true;
    mesh = mesh::Mesh::create(type, iproc, construct_bdrymeshes);
    mesh->readGmsh(ss.str());
    std::cerr << "fini readGmsh " << iproc << "\n";
    mesh->addGeometryObject(meshEnums::MeasureOfCell);
    mesh->addGeometryObject(meshEnums::Normals);
    std::cerr << "fini addGeometryObject " << iproc << "\n";
    std::stringstream ss2;
    ss2 << "tata_" << std::setfill('0') << std::setw(6) << iproc << ".h5";
    mesh->saveH5(ss2.str());
  }
  if(mpirun)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    if(iproc)
    {
      mesh->exchangeInterfaces(nproc);
    }
  }
  if(mpirun) {MPI_Barrier(MPI_COMM_WORLD);}
  if(iproc)
  {
    std::stringstream ss;
    ss << "toto_" << std::setfill('0') << std::setw(6) << iproc << ".h5";
    mesh->saveH5(ss.str());
  }
  if(mpirun) {MPI_Barrier(MPI_COMM_WORLD);}
  if(iproc==0)
  {
    std::cerr << "fini\n";
  }
  // std::shared_ptr<mesh::MeshInterface> mesh = std::make_shared<mesh::TriangleMesh>();
  // mesh->construct_bdrymeshes();
  // mesh->readGmsh("test.msh");
  // mesh->getPlainMesh()->addGeometryObject(meshEnums::MeasureOfCell);
  // mesh->getPlainMesh()->addGeometryObject(meshEnums::Normals);
  // mesh->saveH5("toto.h5");
  // std::cerr << mesh->getInfo() << "\n";
  // mesh->loadH5("mesh.h5");

  if(mpirun)
  {
    std::cerr << "MPI_Finalize " << iproc << "\n";
    MPI_Finalize();
  }
  return 0;
}

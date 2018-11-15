#include  "Alat/sparsitypattern.hpp"
#include  "Mesh/interfacemeshunitsmap.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
InterfaceMeshUnitsMap::~InterfaceMeshUnitsMap() {}
InterfaceMeshUnitsMap::InterfaceMeshUnitsMap(): alat::Map<int,std::shared_ptr<mesh::InterfaceMeshUnit> >(){}
InterfaceMeshUnitsMap::InterfaceMeshUnitsMap( const InterfaceMeshUnitsMap& interfacemeshunitsmap): alat::Map<int,std::shared_ptr<mesh::InterfaceMeshUnit> >(interfacemeshunitsmap)
{
  assert(0);
}
InterfaceMeshUnitsMap& InterfaceMeshUnitsMap::operator=( const InterfaceMeshUnitsMap& interfacemeshunitsmap)
{
  assert(0);
  alat::Map<int,std::shared_ptr<mesh::InterfaceMeshUnit> >::operator=(interfacemeshunitsmap);
  return *this;
}
std::string InterfaceMeshUnitsMap::getClassName() const
{
  return "InterfaceMeshUnitsMap";
}
/*--------------------------------------------------------------------------*/
void InterfaceMeshUnitsMap::exchangeInterfaces(int id, int nproc)
{
  alat::IntSet keys = this->keys();
  int tag = 1;
  MPI_Status status;
  MPI_Request send_request,recv_request;
  double inittime = MPI_Wtime();
  for(alat::IntSet::const_iterator p=keys.begin(); p!=keys.end(); p++)
  {
    if(*p>0) continue;
    int neighbor = -*p;
    // std::cerr << "I want to send from " << id << " to " << neighbor << "\n";
    alat::SparsityPattern sizes_send = (*this)[*p]->getSizes();
    int colsize_send = sizes_send.col().size();
    int rowstartsize_send = sizes_send.rowstart().size();
    int colsize_recieve, rowstartsize_recieve;

    MPI_Isend(&colsize_send, 1, MPI_INT, neighbor, tag, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(&colsize_recieve, 1, MPI_INT, neighbor, tag, MPI_COMM_WORLD, &recv_request);
    MPI_Wait(&recv_request, &status);

    MPI_Isend(&rowstartsize_send, 1, MPI_INT, neighbor, tag, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(&rowstartsize_recieve, 1, MPI_INT, neighbor, tag, MPI_COMM_WORLD, &recv_request);
    MPI_Wait(&recv_request, &status);

    alat::SparsityPattern sizes_receive(colsize_recieve, rowstartsize_recieve);
    MPI_Isend(sizes_send.col().begin(), sizes_send.col().size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(sizes_receive.col().begin(), sizes_receive.col().size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &recv_request);
    MPI_Wait(&recv_request, &status);
    MPI_Isend(sizes_send.rowstart().begin(), sizes_send.rowstart().size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &send_request);
    MPI_Irecv(sizes_receive.rowstart().begin(), sizes_receive.rowstart().size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &recv_request);
    MPI_Wait(&recv_request, &status);
    (*this)[neighbor]->setSizes(sizes_receive);
  }
  for(alat::IntSet::const_iterator p=keys.begin(); p!=keys.end(); p++)
  {
    if(*p>0) continue;
    int neighbor = -*p;
    (*this)[*p]->sendRecv((*this)[neighbor], neighbor);
  }

  double totaltime = MPI_Wtime() - inittime;
  std::cerr << "totaltime = " << totaltime << "\n";
}

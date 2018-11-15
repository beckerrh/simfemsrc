#include  <string>
#include  <cassert>
#include  <iostream>
#include  "timemeshwrapper.hpp"

/*---------------------------------------------------------------------------*/
TimeMeshWrapper::TimeMeshWrapper() : mesh::TimeMesh(){}
TimeMeshWrapper::TimeMeshWrapper(int n, double first, double last) : mesh::TimeMesh(mesh::TimeMeshData(n, first, last)){}
void TimeMeshWrapper::setData(int n, double first, double last)
{
  mesh::TimeMesh::setData(mesh::TimeMeshData(n, first, last));
}
std::ostream &operator<<(std::ostream &os, const TimeMeshWrapper& timemesh)
{
  os << timemesh.getInfo() << "\n";
  return os;
}

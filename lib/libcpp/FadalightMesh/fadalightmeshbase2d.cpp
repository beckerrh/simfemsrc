#include  "FadalightMesh/fadalightmeshbase2d.hpp"
#include  <fstream>

using namespace FadalightMesh;
using namespace std;

/*---------------------------------------------------------*/

template<int NODESPERCELL>
FadalightMeshBase2d<NODESPERCELL>::FadalightMeshBase2d() : FadalightMeshBase<2, NODESPERCELL, NODESPERCELL, 2>( )
{}

template<int NODESPERCELL>
std::string FadalightMeshBase2d<NODESPERCELL>::getClassName() const
{
  return "FadalightMesh::FadalightMeshBase2d";
}

/*---------------------------------------------------------*/

template<int NODESPERCELL>
int FadalightMeshBase2d<NODESPERCELL>::getNEdgesPerCell(int i) const
{
  return 0;
}

/*---------------------------------------------------------*/

template<int NODESPERCELL>
int FadalightMeshBase2d<NODESPERCELL>::getNEdgesPerSide(int i) const
{
  return 0;
}

/*---------------------------------------------------------*/

template<int NODESPERCELL>
int FadalightMeshBase2d<NODESPERCELL>::getVtkType() const
{
  if(NODESPERCELL == 3)
  {
    return 5;
  }
  else
  {
    return 9;
  }
}

/*---------------------------------------------------------*/

template<int NODESPERCELL>
int FadalightMeshBase2d<NODESPERCELL>::getBoundaryVtkType() const
{
  return 3;
}

/*---------------------------------------------------------*/

// triangle mesh
template class FadalightMesh::FadalightMeshBase2d<3>;

// quadrilateral mesh
template class FadalightMesh::FadalightMeshBase2d<4>;

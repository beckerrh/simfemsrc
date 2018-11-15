#include  "FadalightMesh/fadalightmeshbase1d.hpp"
#include  <fstream>

using namespace FadalightMesh;
using namespace std;

/*-------------------------------------------------------*/

template<int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
alat::Node FadalightMeshBase1d<NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodeOfCell(int iK) const
{
  alat::Node v;
  double d = 1.0/ double(NODESPERCELL);
  for(int ii = 0; ii < NODESPERCELL; ii++)
  {
    const alat::Node& vii = getNodeOfCell(iK, ii);
    v.x() += d*vii.x();
  }
  return v;
}

/*---------------------------------------------------------*/
/*---------------------------------------------------------*/

// line mesh
#define NODESPERCELL 2
#define SIDESPERCELL 2
#define NODESPERSIDE 1
template class FadalightMesh::FadalightMeshBase1d<NODESPERCELL, SIDESPERCELL, NODESPERSIDE>;
#undef NODESPERCELL
#undef SIDESPERCELL
#undef NODESPERSIDE

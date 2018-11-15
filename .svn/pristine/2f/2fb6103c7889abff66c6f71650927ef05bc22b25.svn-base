#ifndef __FadalightMesh_FadalightMeshBase2d_h
#define __FadalightMesh_FadalightMeshBase2d_h

#include  "fadalightmeshbase.hpp"

/*---------------------------------------------------------*/

namespace FadalightMesh
{
  template<int NODESPERCELL>
  class FadalightMeshBase2d : public FadalightMeshBase<2, NODESPERCELL, NODESPERCELL, 2>
  {
public:
    FadalightMeshBase2d<NODESPERCELL>( );
    std::string getClassName() const;
    int getNEdgesPerCell(int i) const;
    int getNEdgesPerSide(int i) const;
    int getVtkType() const;
    // std::string getEnsightType() const;
    int getBoundaryVtkType() const;
  };
}

/*---------------------------------------------------------*/

#endif

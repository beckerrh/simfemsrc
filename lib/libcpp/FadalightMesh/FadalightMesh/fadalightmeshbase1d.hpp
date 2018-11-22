#ifndef __FadalightMesh_FadalightMeshBase1d_h
#define __FadalightMesh_FadalightMeshBase1d_h

#include  "fadalightmeshbase.hpp"

/*---------------------------------------------------------*/

namespace FadalightMesh
{
  template<int NODESPERCELL, int SIDESPERCELL, int NODESPERSIDE>
  class FadalightMeshBase1d : public FadalightMeshBase<1, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>
  {
public:

private:

public:

    FadalightMeshBase1d<NODESPERCELL, SIDESPERCELL, NODESPERSIDE>() : FadalightMeshBase<1, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>() {}

    std::string getClassName() const
    {
      return "FadalightMesh::FadalightMeshBase1d";
    }

    const alat::Node& getNode(int i) const
    {
      return FadalightMeshBase<1, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::_nodes[i];
    }

    const alat::Node& getNodeOfSide(int is, int ii) const
    {
      int i =  FadalightMeshBase<1, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodeIdOfSide(is, ii);
      return getNode(i);
    }

    const alat::Node& getNodeOfCell(int iK, int ii) const
    {
      int i =  FadalightMeshBase<1, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNodeIdOfCell(iK, ii);
      assert(i >= 0);
      assert( i < getNNodes() );
      return getNode(i);
    }

    alat::Node getNodeOfCell(int iK) const;

    int getNNodes() const
    {
      return FadalightMeshBase<1, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNNodes();
    }

    int getNCells() const
    {
      return FadalightMeshBase<1, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNCells();
    }

    int getNSides() const
    {
      return FadalightMeshBase<1, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::getNSides();
    }

    // int  getNNodesPerSide() {return NODESPERSIDE;}
  };
}

/*---------------------------------------------------------*/

#endif

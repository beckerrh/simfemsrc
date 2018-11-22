#ifndef  __FadalightMesh_TriangleMesh_h
#define  __FadalightMesh_TriangleMesh_h

#include  "fadalightmeshbase2d.hpp"

/*-----------------------------------------*/

namespace FadalightMesh
{
  class TriangleMesh : public virtual FadalightMesh::FadalightMeshBase2d<3>
  {
public:
    typedef FadalightMesh::FadalightMeshBase2d<3>::Cell Triangle;
    typedef FadalightMesh::FadalightMeshBase2d<3>::Side Edge;
    typedef alat::Node Node;

private:
    mutable Side _S;
    double _ComputeArea(const Triangle& K) const;
    
    //Noeuds de l'element de reference
    alat::Vector<alat::Node> _vhat;
    //Numero des noeuds de chaque face de l'element de reference
    alat::Vector< alat::Vector<int> > _hatnode_id_of_side;
    // alat::Vector< alat::Vector<int> > _hatnode_id_of_edge;
    //Noeuds du patch pour raffinement
    alat::Vector<alat::Node> _patchvhat;
    //Numero des noeuds de chaque face de l'element de volume de reference
    alat::Vector< alat::Vector<int> > _patch_hatnode_id_of_cell;
    // int _nsidespatch, _nedgespatch;

protected:
    const Side& _getSideOfCell(int iK, int ii) const;
    void _hatNodesOfCell();
    void _hatNodeIdOfSide();
    void _hatNodesOfPatch();
    void _patchHatNodeIdOfCell();
    
public:
    ~TriangleMesh();
    TriangleMesh();
    TriangleMesh(const TriangleMesh& trianglemesh);
    TriangleMesh& operator=(const TriangleMesh& trianglemesh);
    std::string getClassName() const;

    int getNodeIdOfSideOfCell(int iK, int iis, int ii) const;

    int getVtkType() const
    {
      return 5;
    }

    void writeMedit(std::string filename) const;
    void ReadTri(std::string filename);
    FadalightMeshEnums::meshtype getType() const;
    std::string getCellType() const;
    int getCouplingOffset(int iS) const;
    void getLocalIndicesOfSidesInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const;
    void getLocalIndicesOfSidesAndDiagonalsInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const;
  };
}

#endif

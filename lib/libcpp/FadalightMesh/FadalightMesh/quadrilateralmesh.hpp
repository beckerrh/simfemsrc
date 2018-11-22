#ifndef __FadalightMesh_QuadrilateralMesh_h
#define __FadalightMesh_QuadrilateralMesh_h

#include  "fadalightmeshbase2d.hpp"


//        2
//     3 -- 2
//   3 |    |  1
//     0 -- 1
//        0

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class QuadrilateralMesh : public FadalightMesh::FadalightMeshBase2d<4>
  {
public:
    typedef FadalightMesh::FadalightMeshBase2d<4>::Cell Quadrilateral;
    typedef FadalightMesh::FadalightMeshBase2d<4>::Side Edge;
    typedef alat::Node Node;

private:
    mutable Side _S;

protected:
    const Side& _getSideOfCell(int iK, int iis) const;

public:
    ~QuadrilateralMesh();
    QuadrilateralMesh();
    QuadrilateralMesh(const QuadrilateralMesh& quadrilateralmesh);
    std::string getClassName() const;
    std::string getCellType() const;

    int getNodeIdOfSideOfCell(int iK, int iis, int ii) const;
    // alat::Vector<QuadrilateralMesh::Quadrilateral>& getQuadrilaterals();
    // const alat::Vector<QuadrilateralMesh::Quadrilateral>& getQuadrilaterals() const;
    // const Quadrilateral& getQuadrilateral(int i) const;
    // const Edge& getEdge(int i) const;
    void read(std::string filename);
    void write(std::string filename) const;
    void writeSimpleMesh(std::string filename, std::string type) const;

    void readQuad(std::string filename);
    void writeQuad(std::string filename) const;
    void check() const;
    // void getMeshSizeForStabilization(double& hs, int iS, int iK, int iil) const;
    int getCouplingOffset(int iS) const;
    FadalightMeshEnums::meshtype getType() const;
    void getLocalIndicesOfSidesInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const;
    void getLocalIndicesOfSidesAndDiagonalsInCell(alat::armaivec& sideindex_a, alat::armaivec& sideindex_e) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

#ifndef __FadalightMesh_HexahedralMesh_h
#define __FadalightMesh_HexahedralMesh_h

#include  "fadalightmeshbase3d.hpp"
#include  "hexahedron.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  //
  //
// #             7-------6
// #             |\      |\
// #             | 4-------5
// #             | |       |
// #             3-|-----2 |
// #              \|      \|
// #               0-------1
//
// #             7-------6                    5
// #             |\   5  |\                4  |
// #             | 4-------5                \ |
// #        2    | |   1   |    3       2-----|--------3
// #             3-|-----2 |                  |\
// #              \|  0   \|                  | \
// #               0-------1                  0  1
//
//

  class HexahedralMesh : public FadalightMesh::FadalightMeshBase3d<8, 6, 12, 4>
  {
public:
    typedef FadalightMesh::FadalightMeshBase3d<8, 6, 12, 4>::Cell Hexahedral;
    typedef FadalightMesh::FadalightMeshBase3d<8, 6, 12, 4>::Side Side;
    typedef alat::Node Node;

    static Hexahedron _hexahedron;

private:
    int getLocalNodeIndiceOfSide(int ii, int isl) const;

private:
    static alat::FixArray<4, int> _ind;
    mutable Edge _E;
    mutable Side _S;
    double _ComputeVolume(const Hexahedral& K) const;

protected:
    const Side& _getSideOfCell(int i, int ii) const;
    const Edge& _getEdgeOfCell(int i, int ii) const;
    const Edge& _getEdgeOfSide(int i, int ii) const;

public:
    HexahedralMesh();
    std::string getClassName() const;
    std::string getCellType() const;
    FadalightMeshEnums::meshtype getType() const;
    alat::Vector<HexahedralMesh::Hexahedral>& getHexahedrals();
    const alat::Vector<HexahedralMesh::Hexahedral>& getHexahedrals() const;
    const Hexahedral& getHexahedral(int i) const;
    void read(std::string filename);
    void write(std::string filename) const;
    void readHex(std::string filename);
    void writeHex(std::string filename) const;

     void writeVtk(std::string filename) const;
    // void writeBoundaryVtk(std::string filename) const;
    Side getSideOfCell(int i, int ii) const;
    int getNodeIdOfSideOfCell(int iK, int iis, int ii) const;
    // void getMeshSizeForStabilization(double& hs, int iS, int iK, int iil) const;
    int getCouplingOffset(int iS) const;
    int getBoundaryVtkType() const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

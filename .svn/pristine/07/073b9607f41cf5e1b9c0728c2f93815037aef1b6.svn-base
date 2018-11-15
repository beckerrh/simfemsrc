#ifndef __FadalightMesh_FadalightMeshBase3d_h
#define __FadalightMesh_FadalightMeshBase3d_h

#include  "fadalightmeshbase.hpp"

/*---------------------------------------------------------*/

namespace FadalightMesh
{
  template<int NODESPERCELL, int SIDESPERCELL, int EDGESPERCELL, int NODESPERSIDE>
  class FadalightMeshBase3d : public FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>
  {
public:
    typedef CellBase<2> Edge;
    typedef CellBase<EDGESPERCELL> EdgeCell;
    typedef CellBase<NODESPERSIDE> EdgeSide;
    typedef CellBase<SIDESPERCELL> SideCell;
    typedef CellBase<NODESPERSIDE> Side;
    // typedef std::map<Side, int> BoundarySideToColor;
    typedef typename FadalightMeshBase<3, NODESPERCELL, SIDESPERCELL, NODESPERSIDE>::BoundarySideToColor BoundarySideToColor;

private:
    static Edge _dummyedge;
    alat::Vector<Edge>   _edges;
    alat::Vector<EdgeCell>   _edges_of_cells;
    alat::Vector<EdgeSide>   _edges_of_sides;

protected:
    virtual const Edge& _getEdgeOfCell(int i, int ii) const;
    virtual const Edge& _getEdgeOfSide(int i, int ii) const;

public:
    FadalightMeshBase3d<NODESPERCELL, SIDESPERCELL, EDGESPERCELL, NODESPERSIDE>( );
    alat::Vector<Edge>& getEdges();
    alat::Vector<EdgeCell>&  getEdgesOfCells();
    alat::Vector<EdgeSide>&  getEdgesOfSides();
    const alat::Vector<Edge>&   getEdges() const;
    const alat::Vector<EdgeCell>&  getEdgesOfCells() const;
    const alat::Vector<EdgeSide>&  getEdgesOfSides() const;
    const Edge&   getEdge(int i) const;
    const EdgeCell&  getEdgesOfCell(int i) const;
    const EdgeSide&  getEdgesOfSides(int i) const;
    std::string getClassName() const;
    int getNEdgesPerCell(int i) const;
    int getNEdgesPerSide(int i) const;
    int getNEdges() const;
    int getNodeIdOfEdge(int i, int ii) const;
    int getEdgeIdOfCell(int i, int ii) const;
    int getEdgeIdOfSide(int i, int ii) const;
    alat::Node getNodeOfEdge(int ie) const;
    void readFadalightMesh(const std::string& basefilename);
    void writeFadalightMesh(const std::string& basefilename, arma::file_type datatype = arma::arma_binary) const;
    void constructSidesFromCells(BoundarySideToColor& bsides, int color_default = 0);
    // std::string getEnsightType() const;
    int getVtkType() const;
  };
}

/*---------------------------------------------------------*/

#endif

#ifndef __Mesh_MeshVisitorTriangle_hpp
#define __Mesh_MeshVisitorTriangle_hpp

#include  "meshvisitorinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class MeshVisitorTriangle : public MeshVisitorInterface
  {
  public:
    ~MeshVisitorTriangle();
    MeshVisitorTriangle();
    MeshVisitorTriangle( const MeshVisitorTriangle& meshvisitortriangle);
    MeshVisitorTriangle& operator=( const MeshVisitorTriangle& meshvisitortriangle);
    std::string getClassName() const;
    meshEnums::meshtype getType() const;
    std::unique_ptr<MeshVisitorInterface> clone() const;

    std::unique_ptr<MeshVisitorInterface> newBoundaryVisitor() const;
    void set_dimensions(int& ndim, int& nnodes_per_cell, int& nedges_per_cell, int& nsides_per_cell, int& nnodes_per_side) const;
    double computeMeasureOfCell(int iK, const mesh::MeshUnitInterface* mesh) const;
    void computeNormal(MeshUnitInterface::Normal normal, const mesh::MeshUnitInterface* mesh, const alat::Node& xS, const MeshUnitInterface::Side S, int iK0, int ii0, int iK1, int ii1) const;

  };
}

/*--------------------------------------------------------------------------*/
#endif

#ifndef __Mesh_MeshVisitorPoint_hpp
#define __Mesh_MeshVisitorPoint_hpp

#include  "meshvisitorinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class MeshVisitorPoint : public MeshVisitorInterface
  {
  public:
    ~MeshVisitorPoint();
    MeshVisitorPoint();
    MeshVisitorPoint( const MeshVisitorPoint& meshvisitorpoint);
    MeshVisitorPoint& operator=( const MeshVisitorPoint& meshvisitorpoint);
    std::string getClassName() const;
    meshEnums::meshtype getType() const;
    std::unique_ptr<MeshVisitorInterface> clone() const;

    std::unique_ptr<MeshVisitorInterface> newBoundaryVisitor() const;
    void set_dimensions(int& ndim, int& nnodes_per_cell, int& nedges_per_cell, int& nsides_per_cell, int& nnodes_per_side) const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

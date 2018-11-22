#ifndef __Mesh_MeshVisitorInterface_hpp
#define __Mesh_MeshVisitorInterface_hpp

#include  "Alat/interfacebase.hpp"
#include  "Mesh/meshunitinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class MeshVisitorInterface : public alat::InterfaceBase
  {
  public:
    ~MeshVisitorInterface();
    MeshVisitorInterface();
    MeshVisitorInterface( const MeshVisitorInterface& meshvisitorinterface);
    MeshVisitorInterface& operator=( const MeshVisitorInterface& meshvisitorinterface);
    std::string getClassName() const;
    virtual meshEnums::meshtype getType() const=0;
    virtual std::unique_ptr<MeshVisitorInterface> clone() const=0;

    virtual std::unique_ptr<MeshVisitorInterface> newBoundaryVisitor() const=0;
    virtual void set_dimensions(int& ndim, int& nnodes_per_cell, int& nedges_per_cell, int& nsides_per_cell, int& nnodes_per_side) const=0;
    virtual double computeMeasureOfCell(int iK, const mesh::MeshUnitInterface* mesh) const;
    virtual void computeNormal(MeshUnitInterface::Normal normal, const mesh::MeshUnitInterface* mesh, const alat::Node& xS, const MeshUnitInterface::Side S, int iK0, int ii0, int iK1, int ii1) const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

#ifndef __Mesh_BoundaryMeshUnit_hpp
#define __Mesh_BoundaryMeshUnit_hpp

#include  "meshunit.hpp"
#include  "boundarymeshinfo.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class BoundaryMeshUnit : public virtual MeshUnit
  {
  protected:
    BoundaryMeshInfo _boundarymeshinfo;
    const MeshUnitInterface* _parentmesh;

  public:
    ~BoundaryMeshUnit();
    BoundaryMeshUnit();
    BoundaryMeshUnit( const BoundaryMeshUnit& boundarymeshunit);
    void init(std::shared_ptr<const mesh::MeshVisitorInterface> visitor, const MeshUnitInterface* parentmesh);
    BoundaryMeshUnit& operator=( const BoundaryMeshUnit& boundarymeshunit);
    std::string getClassName() const;
    const MeshUnitInterface* getParentMesh() const;
    const BoundaryMeshInfo& getBoundaryMeshInfo() const;
    BoundaryMeshInfo& getBoundaryMeshInfo();
    void saveH5(const arma::hdf5_name& spec) const;
    void loadH5(const arma::hdf5_name& spec);
  };
}

/*--------------------------------------------------------------------------*/
#endif

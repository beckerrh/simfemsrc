#ifndef __mesh_InterfaceMeshUnit_hpp
#define __mesh_InterfaceMeshUnit_hpp

#include  "meshunit.hpp"
#include  "interfacemeshinfo.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class InterfaceMeshUnit : public MeshUnit
  {
  protected:
    InterfaceMeshInfo _interfacemeshinfo;
    const MeshUnitInterface* _plainmesh;

  public:
    ~InterfaceMeshUnit();
    InterfaceMeshUnit();
    InterfaceMeshUnit( const InterfaceMeshUnit& interfacemeshunit);
    void init(std::shared_ptr<const mesh::MeshVisitorInterface> visitor, const MeshUnitInterface* plainmesh);
    InterfaceMeshUnit& operator=( const InterfaceMeshUnit& interfacemeshunit);
    std::string getClassName() const;
    const MeshUnitInterface* getPlainMesh() const;
    const InterfaceMeshInfo& getInterfaceMeshInfo() const;
    InterfaceMeshInfo& getInterfaceMeshInfo();
    void saveH5(const arma::hdf5_name& spec) const;
    void loadH5(const arma::hdf5_name& spec);
  };
}

/*--------------------------------------------------------------------------*/
#endif

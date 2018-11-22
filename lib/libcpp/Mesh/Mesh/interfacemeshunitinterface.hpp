#ifndef __Mesh_InterfaceMeshUnitInterface_hpp
#define __Mesh_InterfaceMeshUnitInterface_hpp

#include  "meshunitinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class InterfaceMeshUnitInterface : public MeshUnitInterface
  {
  private:
  protected:
  public:
    ~InterfaceMeshUnitInterface();
    InterfaceMeshUnitInterface();
    InterfaceMeshUnitInterface( const InterfaceMeshUnitInterface& interfacemeshunitinterface);
    InterfaceMeshUnitInterface& operator=( const InterfaceMeshUnitInterface& interfacemeshunitinterface);
    std::string getClassName() const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

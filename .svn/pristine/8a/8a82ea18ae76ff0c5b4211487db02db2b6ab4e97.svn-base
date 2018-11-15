#ifndef __Mesh_InterfaceMeshUnitsMap_hpp
#define __Mesh_InterfaceMeshUnitsMap_hpp

#include  "Alat/map.hpp"
#include  "Mesh/interfacemeshunit.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class InterfaceMeshUnitsMap : public alat::Map<int,std::shared_ptr<mesh::InterfaceMeshUnit> >
  {
  private:
  protected:
  public:
    ~InterfaceMeshUnitsMap();
    InterfaceMeshUnitsMap();
    InterfaceMeshUnitsMap( const InterfaceMeshUnitsMap& interfacemeshunitsmap);
    InterfaceMeshUnitsMap& operator=( const InterfaceMeshUnitsMap& interfacemeshunitsmap);
    std::string getClassName() const;

    void exchangeInterfaces(int id, int nproc);
  };
}

/*--------------------------------------------------------------------------*/
#endif

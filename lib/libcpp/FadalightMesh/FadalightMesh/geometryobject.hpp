#ifndef __FadalightMesh_GeometryObject_h
#define __FadalightMesh_GeometryObject_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/armadillo.hpp"
// #include  "FadalightMesh/meshinterface.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class MeshInterface;

  class GeometryObject : public virtual alat::InterfaceBase
  {
protected:
    std::string getInterfaceName() const;

public:
    ~GeometryObject();
    GeometryObject();
    GeometryObject( const GeometryObject& geometryobject);
    GeometryObject& operator=( const GeometryObject& geometryobject);

    virtual void load(std::string filename) = 0;
    virtual void save(std::string filename, arma::file_type datatype = arma::arma_binary) const = 0;
    virtual void constructGeometryObject(const FadalightMesh::MeshInterface* mesh);
  };
}

/*--------------------------------------------------------------------------*/

#endif

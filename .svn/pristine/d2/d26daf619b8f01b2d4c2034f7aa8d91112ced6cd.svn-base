#ifndef __Mesh_InterfaceMeshInfo_hpp
#define __Mesh_InterfaceMeshInfo_hpp

#include  "geometryobject.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class InterfaceMeshInfo : public GeometryObject
  {
  public:
    alat::armaivec  _nodes_to_plain;
    alat::armaivec  _sides_to_plain;
    alat::armaivec  _cells_to_plain;

public:
    ~InterfaceMeshInfo();
    InterfaceMeshInfo();
    InterfaceMeshInfo( const InterfaceMeshInfo& interfacemeshinfo);
    InterfaceMeshInfo& operator=( const InterfaceMeshInfo& interfacemeshinfo);
    std::string getClassName() const;
    std::unique_ptr<GeometryObject> clone() const;

    const alat::armaivec&  getNodesToPlain() const;
    alat::armaivec&  getNodesToPlain();
    const alat::armaivec&  getCellsToPlain() const;
    alat::armaivec&  getCellsToPlain();
    const alat::armaivec&  getSidesToPlain() const;
    alat::armaivec&  getSidesToPlain();

    alat::armaivec getSizes() const;
    void setSizes(alat::armaivec::const_iterator sizes);
    void send(int neighbor, int tag) const;
    void recv(int neighbor, int tag);
    void loadH5(const arma::hdf5_name& spec);
    void saveH5(const arma::hdf5_name& spec) const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

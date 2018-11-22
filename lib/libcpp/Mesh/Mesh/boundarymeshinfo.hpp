#ifndef __Mesh_BoundaryMeshInfo_hpp
#define __Mesh_BoundaryMeshInfo_hpp

#include  "geometryobject.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class BoundaryMeshInfo : public GeometryObject
  {
  public:
    alat::armaivec  _nodes_to_parent;
    alat::armaimat  _cells_to_parent;

public:
    ~BoundaryMeshInfo();
    BoundaryMeshInfo();
    BoundaryMeshInfo( const BoundaryMeshInfo& sidesandcells);
    BoundaryMeshInfo& operator=( const BoundaryMeshInfo& sidesandcells);
    std::string getClassName() const;
    std::unique_ptr<GeometryObject> clone() const;

    const alat::armaivec&  getNodesToParent() const;
    alat::armaivec&  getNodesToParent();
    const alat::armaimat&  getCellsToParent() const;
    alat::armaimat&  getCellsToParent();

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

#ifndef __Mesh_NodesCellsWeight_hpp
#define __Mesh_NodesCellsWeight_hpp

#include  "Alat/armadillo.hpp"
#include  "geometryobject.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class NodesCellsWeight : public GeometryObject
  {
  protected:
    arma::vec _weights;

  public:
    ~NodesCellsWeight();
    NodesCellsWeight();
    NodesCellsWeight( const NodesCellsWeight& nodescellsweight);
    NodesCellsWeight& operator=( const NodesCellsWeight& nodescellsweight);
    std::string getClassName() const;
    std::unique_ptr<GeometryObject> clone() const;

    const arma::vec& get() const;

    alat::armaivec getSizes() const;
    void setSizes(alat::armaivec::const_iterator sizes);
    void send(int neighbor, int tag) const;
    void recv(int neighbor, int tag);
    void loadH5(const arma::hdf5_name& spec);
    void saveH5(const arma::hdf5_name& spec) const;
    void construct(const mesh::MeshUnitInterface* mesh, const GeometryConstructorInterface* geometryconstructor);
  };
}

/*--------------------------------------------------------------------------*/
#endif

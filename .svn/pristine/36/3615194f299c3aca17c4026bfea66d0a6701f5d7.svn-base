#ifndef __Mesh_Normals_h
#define __Mesh_Normals_h

#include  "Alat/vector.hpp"
#include  "Mesh/geometryobject.hpp"

/*--------------------------------------------------------------------------*/

namespace mesh
{
  class Normals : public virtual GeometryObject
  {
protected:
    arma::mat  _normals;
    arma::fmat _sigma;

public:
    ~Normals();
    Normals();
    Normals( const Normals& geometryobject);
    Normals& operator=( const Normals& geometryobject);
    std::string getClassName() const;
    std::unique_ptr<GeometryObject> clone() const;

    const arma::mat& getNormals() const;
    const arma::fmat& getSigma() const;

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

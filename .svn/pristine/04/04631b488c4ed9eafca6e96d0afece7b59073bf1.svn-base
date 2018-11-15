#ifndef __Mesh_GeometryObject_h
#define __Mesh_GeometryObject_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/set.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class SparsityPattern;
}
namespace mesh
{
  class MeshUnitInterface;
}
namespace mesh
{
  class GeometryConstructorInterface : public virtual alat::InterfaceBase{
  public: std::string getInterfaceName() const {return "GeometryConstructorInterface";}
  };
  class TrivialGeometryConstructor : public GeometryConstructorInterface{
    public: std::string getClassName() const {return "TrivialGeometryConstructor";}
  };
  static const TrivialGeometryConstructor trivialgeometryconstructor;

  class GeometryObject : public virtual alat::InterfaceBase
  {
protected:
    std::string getInterfaceName() const;

public:
    ~GeometryObject();
    GeometryObject();
    GeometryObject( const GeometryObject& geometryobject);
    GeometryObject& operator=( const GeometryObject& geometryobject);
    virtual std::unique_ptr<GeometryObject> clone() const=0;

    virtual alat::armaivec getSizes() const=0;
    virtual void setSizes(alat::armaivec::const_iterator sizes)=0;
    virtual void send(int neighbor, int tag) const=0;
    virtual void recv(int neighbor, int tag)=0;
    virtual void loadH5(const arma::hdf5_name& spec) = 0;
    virtual void saveH5(const arma::hdf5_name& spec) const = 0;
    virtual void construct(const mesh::MeshUnitInterface* mesh, const GeometryConstructorInterface* geometryconstructor=&trivialgeometryconstructor);
  };
}

/*--------------------------------------------------------------------------*/

#endif

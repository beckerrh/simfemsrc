#ifndef __Mesh_MeasureOfCell_h
#define __Mesh_MeasureOfCell_h

#include  "Alat/armadillo.hpp"
#include  "Mesh/geometryobject.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class MeasureOfCell : public virtual GeometryObject
  {
  protected:
    arma::vec _measureofcell;
    // arma::mat _centerofcell;

  public:
    ~MeasureOfCell();
    MeasureOfCell();
    MeasureOfCell( const MeasureOfCell& geometryobject);
    MeasureOfCell& operator=( const MeasureOfCell& geometryobject);
    std::string getClassName() const;
    std::unique_ptr<GeometryObject> clone() const;

    const arma::vec& getMeasureOfCell() const;
    // const arma::mat& getCenterOfCell() const;

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

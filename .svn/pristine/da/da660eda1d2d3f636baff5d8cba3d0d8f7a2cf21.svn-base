#ifndef __Mesh_BoundaryInformation_hpp
#define __Mesh_BoundaryInformation_hpp

#include  "geometryobject.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class BoundaryInformation : public GeometryObject
  {
  protected:
    alat::armaimat _cells_on_bdry_of_plain;

  public:
    ~BoundaryInformation();
    BoundaryInformation();
    BoundaryInformation( const BoundaryInformation& boundarymeshinformation);
    BoundaryInformation& operator=( const BoundaryInformation& boundarymeshinformation);
    std::string getClassName() const;
    std::unique_ptr<GeometryObject> clone() const;

    const alat::armaimat& getCellsOnBdryOfPlain() const;
    alat::armaimat& getCellsOnBdryOfPlain();

    void set_size(int n);
    int size() const;
    int cellid(int i) const;
    int sideid(int i) const;
    int sidelocalid(int i) const;

    alat::armaivec getSizes() const;
    void setSizes(alat::armaivec::const_iterator sizes);
    void send(int neighbor, int tag) const;
    void recv(int neighbor, int tag);
    void loadH5(const arma::hdf5_name& spec);
    void saveH5(const arma::hdf5_name& spec) const;
  };
  std::ostream& operator<<(std::ostream& os, const BoundaryInformation& boundaryinformation);
}

/*--------------------------------------------------------------------------*/
#endif

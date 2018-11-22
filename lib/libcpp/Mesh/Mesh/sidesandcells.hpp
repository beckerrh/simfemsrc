#ifndef __Mesh_SidesAndCells_hpp
#define __Mesh_SidesAndCells_hpp

#include  "geometryobject.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class SidesAndCells : public GeometryObject
  {
  public:
    alat::armaimat  _nodes_of_sides;
    alat::armaimat  _sides_of_cells;
    alat::armaimat  _cells_of_sides;

  public:
    ~SidesAndCells();
    SidesAndCells();
    SidesAndCells( const SidesAndCells& sidesandcells);
    SidesAndCells& operator=( const SidesAndCells& sidesandcells);
    std::string getClassName() const;
    std::unique_ptr<GeometryObject> clone() const;

    const alat::armaimat& getNodesOfSides() const;
    alat::armaimat& getNodesOfSides();
    const alat::armaimat& getSidesOfCells() const;
    alat::armaimat& getSidesOfCells();
    const alat::armaimat& getCellsOfSides() const;
    alat::armaimat& getCellsOfSides();

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

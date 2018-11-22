#ifndef __Mesh_EdgesAndCells_hpp
#define __Mesh_EdgesAndCells_hpp

#include  "geometryobject.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class EdgesAndCells : public GeometryObject
  {
  public:
    alat::armaimat  _nodes_of_edges;
    alat::armaimat  _edges_of_cells;
    alat::armaicube  _localnodes_of_edges_in_cells;

  public:
    ~EdgesAndCells();
    EdgesAndCells();
    EdgesAndCells( const EdgesAndCells& edgesandcells);
    EdgesAndCells& operator=( const EdgesAndCells& edgesandcells);
    std::string getClassName() const;
    std::unique_ptr<GeometryObject> clone() const;

    const alat::armaimat& getNodesOfEdges() const;
    alat::armaimat& getNodesOfEdges();
    const alat::armaimat& getEdgesOfCells() const;
    alat::armaimat& getEdgesOfCells();
    const alat::armaicube& getLocalnodesOfEdgesOfCells() const;
    alat::armaicube& getLocalnodesOfEdgesOfCells();

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

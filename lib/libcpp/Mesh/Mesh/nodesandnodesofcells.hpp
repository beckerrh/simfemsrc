#ifndef __Mesh_NodesAndNodesOfCells_hpp
#define __Mesh_NodesAndNodesOfCells_hpp

#include  "geometryobject.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class NodesAndNodesOfCells : public GeometryObject
  {
  public:
    arma::mat  _nodes;
    alat::armaimat  _nodes_of_cells;

  public:
    ~NodesAndNodesOfCells();
    NodesAndNodesOfCells();
    NodesAndNodesOfCells( const NodesAndNodesOfCells& nodesandnodesofcells);
    NodesAndNodesOfCells& operator=( const NodesAndNodesOfCells& nodesandnodesofcells);
    std::string getClassName() const;
    std::unique_ptr<GeometryObject> clone() const;

    const arma::mat& getNodes() const;
    arma::mat& getNodes();

    const alat::armaimat& getNodesOfCells() const;
    alat::armaimat& getNodesOfCells();

    alat::armaivec getSizes() const;
    void setSizes(alat::armaivec::const_iterator sizes);
    void send(int neighbor, int tag) const;
    void recv(int neighbor, int tag);
    // void sendRecv(NodesAndNodesOfCells& recieve, int neighbor, int tag) const;
    void loadH5(const arma::hdf5_name& spec);
    void saveH5(const arma::hdf5_name& spec) const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

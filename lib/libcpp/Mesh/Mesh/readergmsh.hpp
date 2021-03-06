#ifndef __Mesh_ReaderGmsh_h
#define __Mesh_ReaderGmsh_h

#include  "Alat/armadillo.hpp"
#include  "Alat/interfacebase.hpp"
#include  "Alat/sparsitypattern.hpp"
#include  "Mesh/meshunitinterface.hpp"

/*---------------------------------------------------------*/
namespace mesh
{
  class ReaderGmsh : public virtual alat::InterfaceBase
  {
protected:
    // std::string filename;
    void _read_nodes_binaray(std::ifstream& file);
    void _read_nodes_ascii(std::ifstream& file);
    void _read_elements_binaray(std::ifstream& file);
    void _read_elements_ascii(std::ifstream& file);

    void _setMeshAndInterfaces(MeshUnitInterface* mesh);

    bool _binary;
    int _npartsmax;
    alat::armaivec _size_of_elem;

public:
    alat::Map<int, alat::armaimat> color_to_bdry1, color_to_bdry2, color_to_bdry3;
    alat::Map<int, alat::armaivec> partition_to_cells;
    arma::mat _nodes;
    alat::SparsityPattern _elems;
    alat::IntMap _nodeid2id, _cellid2id;


    ~ReaderGmsh();
    ReaderGmsh();
    std::string getClassName() const;

    void read(std::string filename);
    void setMesh(MeshUnitInterface* mesh);
  };
}

/*---------------------------------------------------------*/

#endif

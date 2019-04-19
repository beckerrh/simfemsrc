#ifndef __Mesh_ReaderGmshNew_h
#define __Mesh_ReaderGmshNew_h

#include  "Alat/armadillo.hpp"
#include  "Alat/interfacebase.hpp"
#include  "Alat/sparsitypattern.hpp"
#include  "Alat/vector.hpp"
#include  "Mesh/meshunitinterface.hpp"

/*---------------------------------------------------------*/
namespace mesh
{
  struct EntityTagPoint {
    alat::armaivec tag;
    alat::armamat coord;
    alat::Vector<alat::armaivec> physicaltag;
  };
  struct EntityTag {
    alat::armaivec tag;
    alat::armamat coordmin, coordmax;
    alat::Vector<alat::armaivec> physicaltag;
    alat::Vector<alat::armaivec> physicaltaglower;
  };
  class ReaderGmshNew : public virtual alat::InterfaceBase
  {
protected:
    // std::string filename;
    void _read_nodes_binaray(std::ifstream& file);
    void _read_nodes_ascii(std::ifstream& file);
    void _read_elements_binaray(std::ifstream& file);
    void _read_elements_ascii(std::ifstream& file);

    void _read_physicalnames_binaray(std::ifstream& file);
    void _read_physicalnames_ascii(std::ifstream& file);
    void _read_entities_binaray(std::ifstream& file);
    void _read_entities_ascii(std::ifstream& file);
    void _read_partitionedentities_binaray(std::ifstream& file);
    void _read_partitionedentities_ascii(std::ifstream& file);
    void _read_periodic_binaray(std::ifstream& file);
    void _read_periodic_ascii(std::ifstream& file);
    void _read_ghostelements_binaray(std::ifstream& file);
    void _read_ghostelements_ascii(std::ifstream& file);
    void _read_nodedata_binaray(std::ifstream& file);
    void _read_nodedata_ascii(std::ifstream& file);
    void _read_elementdata_binaray(std::ifstream& file);
    void _read_elementdata_ascii(std::ifstream& file);
    void _read_elementnodedata_binaray(std::ifstream& file);
    void _read_elementnodedata_ascii(std::ifstream& file);
    void _read_interpolationscheme_binaray(std::ifstream& file);
    void _read_interpolationscheme_ascii(std::ifstream& file);

    void _setMeshAndInterfaces(MeshUnitInterface* mesh);

    bool _binary;
    int _npartsmax;
    alat::armaivec _size_of_elem;

    EntityTagPoint _pointtag;
    EntityTag _curvetag, _surfacetag, _volumetag;


public:
    alat::Map<int, alat::armaimat> color_to_bdry1, color_to_bdry2, color_to_bdry3;
    alat::Map<int, alat::armaivec> partition_to_cells;
    arma::mat _nodes;
    alat::SparsityPattern _elems;
    alat::IntMap _nodeid2id, _cellid2id;


    ~ReaderGmshNew();
    ReaderGmshNew();
    std::string getClassName() const;

    void read(std::string filename);
    void setMesh(MeshUnitInterface* mesh);
  };
}

/*---------------------------------------------------------*/

#endif

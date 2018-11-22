
#ifndef __Mesh_EdgesConstructor_h
#define __Mesh_EdgesConstructor_h

#include  "Alat/fixarray.hpp"
#include  "Alat/map.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  "Alat/armadillo.hpp"

/*---------------------------------------------------------*/
namespace mesh
{
  class EdgesConstructor
  {
  protected:
    typedef alat::FixArray<2, int> FixEdge;
    typedef alat::Map<FixEdge, alat::armaimat >  EdgeToInfo;

    MeshUnitInterface* _mesh;

  public:
    ~EdgesConstructor();
    EdgesConstructor(MeshUnitInterface* mesh);
    void constructEdgesFromCells(bool debug=false);
  };
}

/*---------------------------------------------------------*/

#endif

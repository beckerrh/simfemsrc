
#ifndef __Mesh_SidesConstructor_h
#define __Mesh_SidesConstructor_h

#include  "Alat/fixarray.hpp"
#include  "Alat/map.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  "Alat/armadillo.hpp"

/*---------------------------------------------------------*/
namespace mesh
{
  template<int NODESPERSIDE>
  class SidesConstructor
  {
  protected:
    typedef alat::FixArray<NODESPERSIDE, int> FixSide;
    typedef alat::Map<FixSide, int> BoundarySideToColor;
    typedef alat::Map<FixSide, alat::FixArray<2, int> >  SideToInfo;

    MeshUnitInterface* _mesh;

  public:
    ~SidesConstructor<NODESPERSIDE>();
    SidesConstructor<NODESPERSIDE>(MeshUnitInterface* mesh);
    void constructSidesFromCells(const MeshUnitInterface::BoundarySideToColorForSidesConstructor& bsides,
      const alat::Map<int, alat::armaivec>* partition_to_cells=NULL, int color_default = 0, bool debug=false);
  };
}

/*---------------------------------------------------------*/

#endif

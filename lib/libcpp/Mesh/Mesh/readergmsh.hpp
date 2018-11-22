
#ifndef __Mesh_ReaderGmsh_h
#define __Mesh_ReaderGmsh_h

#include  "Alat/interfacebase.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  "Alat/armadillo.hpp"

/*---------------------------------------------------------*/
namespace mesh
{
  class ReaderGmsh : public virtual alat::InterfaceBase
  {
  protected:
    MeshUnitInterface* _mesh;

  public:
    alat::Map<int, alat::armaimat> color_to_bdry1, color_to_bdry2, color_to_bdry3;

    ~ReaderGmsh();
    ReaderGmsh(MeshUnitInterface* mesh);
    std::string getClassName() const;

    void read(std::string filename);
  };
}

/*---------------------------------------------------------*/

#endif

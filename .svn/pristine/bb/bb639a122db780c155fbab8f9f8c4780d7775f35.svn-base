#ifndef  __MeshWrap_MeshWrapper_h
#define  __MeshWrap_MeshWrapper_h

#include  "Mesh/mesh.hpp"

/*---------------------------------------------------------------------------*/
class MeshWrapper : public mesh::Mesh
{
public:
  MeshWrapper(const mesh::Mesh& mesh) : mesh::Mesh(mesh){}
};

/*---------------------------------------------------------------------------*/
std::shared_ptr<MeshWrapper> create(std::string const& s, int partion_id, bool construct_bdrymeshes);
// std::shared_ptr<mesh::Mesh> create(std::string const& s, int partion_id, bool construct_bdrymeshes);

void wrapMesh();



#endif

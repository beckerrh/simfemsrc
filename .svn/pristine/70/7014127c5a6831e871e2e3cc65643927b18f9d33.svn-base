#include  <string>
#include  <cassert>
#include  <iostream>
#include  "meshwrapper.hpp"

/*---------------------------------------------------------------------------*/
std::shared_ptr<MeshWrapper> create(std::string const& s, int partion_id, bool construct_bdrymeshes)
{
  std::shared_ptr<mesh::MeshInterface> meshi = mesh::Mesh::create(s, partion_id, construct_bdrymeshes);
  std::shared_ptr<mesh::Mesh> meshptr = std::dynamic_pointer_cast<mesh::Mesh>(meshi);
  return std::make_shared<MeshWrapper>(*meshptr);
}
//
// /*---------------------------------------------------------------------------*/
// std::shared_ptr<mesh::Mesh> create(std::string const& s, int partion_id, bool construct_bdrymeshes)
// {
//   std::shared_ptr<mesh::MeshInterface> meshi = mesh::Mesh::create(s, partion_id, construct_bdrymeshes);
//   std::shared_ptr<mesh::Mesh> meshptr = std::dynamic_pointer_cast<mesh::Mesh>(meshi);
//   return meshptr;
// }

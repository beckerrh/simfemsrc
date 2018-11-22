#include  "Solvers/modelinterface.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
ModelInterface::~ModelInterface() {}
ModelInterface::ModelInterface(): alat::InterfaceBase(){}
ModelInterface::ModelInterface( const ModelInterface& modelinterface): alat::InterfaceBase(modelinterface)
{
assert(0);
}
ModelInterface& ModelInterface::operator=( const ModelInterface& modelinterface)
{
  assert(0);
  alat::InterfaceBase::operator=(modelinterface);
  return *this;
}
std::string ModelInterface::getClassName() const
{
  return "ModelInterface";
}
/*--------------------------------------------------------------------------*/
void ModelInterface::initModel(const mesh::MeshUnitInterface* mesh, const solvers::Parameters& parameters)
{
  _mesh=mesh;
}

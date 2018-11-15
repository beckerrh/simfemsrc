#include  "Solvers/model.hpp"
#include  "Solvers/feminterface.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
Model::~Model() {}
Model::Model(): ModelInterface(){}
Model::Model( const Model& model): ModelInterface(model)
{
assert(0);
}
Model& Model::operator=( const Model& model)
{
  assert(0);
  ModelInterface::operator=(model);
  return *this;
}
std::string Model::getClassName() const
{
  return "Model";
}
std::string Model::getInfo() const
{
  return "";
}

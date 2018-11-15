#include  "model.hpp"
#include  "Solvers/solverinterface.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
Model::~Model() {}
Model::Model(): solvers::ModelInterface(){}
Model::Model( const Model& model): solvers::ModelInterface(model)
{
assert(0);
}
Model& Model::operator=( const Model& model)
{
  assert(0);
  solvers::ModelInterface::operator=(model);
  return *this;
}
std::string Model::getClassName() const
{
  return "Model";
}
std::string Model::getInfo() const
{
  std::stringstream ss;
  ss << "\t a="<<a<<" b="<<b<<" k0="<<k0<<" k1="<<k1<<"\n";
  return ss.str();
}

/*--------------------------------------------------------------------------*/
void Model::initModel(const mesh::MeshUnitInterface* mesh, const solvers::Parameters& parameters)
{
  a = parameters.doubles["a"];
  b = parameters.doubles["b"];
  k0 = parameters.doubles["k0"];
  k1 = parameters.doubles["k1"];
}

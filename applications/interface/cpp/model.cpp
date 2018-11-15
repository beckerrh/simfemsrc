#include  "Solvers/model.hpp"
#include  "Solvers/feminterface.hpp"
#include  "Solvers/solverinterface.hpp"
#include  "model.hpp"
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
  ss << "\t _kin="<<_kin<<" _kex="<<_kex<<" _xgamma="<<_xgamma;
  return ss.str();
}

/*--------------------------------------------------------------------------*/
void Model::initModel(const mesh::MeshUnitInterface* mesh, const solvers::Parameters& parameters)
{
  solvers::ModelInterface::initModel(mesh, parameters);
  assert(mesh);
  _kin = parameters.doubles["kin"];
  _kex = parameters.doubles["kex"];
  _xgamma = parameters.doubles["xgamma"];
}

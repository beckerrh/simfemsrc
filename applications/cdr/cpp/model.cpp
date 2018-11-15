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
  ss << "\t alpha="<<_alpha<<" diff="<<_diff<<" beta="<<_beta->getClassName()<<"\n";
  return ss.str();
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
class BetaZero : public solvers::InitialConditionInterface
{
public:
  std::string getClassName()const {return "BetaZero";}
  void operator()(arma::vec& beta, double x, double y, double z, double t)const{
    for(int i=0;i<_dim;i++) {beta[i] = 0.0;}
  }
};
class BetaConstant : public solvers::InitialConditionInterface
{
public:
  std::string getClassName()const {return "BetaConstant";}
  void operator()(arma::vec& beta, double x, double y, double z, double t)const{
    for(int i=0;i<_dim;i++) {beta[i] = (double) _dim;}
  }
};
class BetaEast : public solvers::InitialConditionInterface
{
public:
  std::string getClassName()const {return "BetaEast";}
  void operator()(arma::vec& beta, double x, double y, double z, double t)const{
    assert(beta.size()==3);
    beta[0]=1.0;
    beta[1]=0.0;
    beta[2]=0.0;
  }
};

/*--------------------------------------------------------------------------*/
void Model::initModel(const mesh::MeshUnitInterface* mesh, const solvers::Parameters& parameters)
{
  solvers::ModelInterface::initModel(mesh, parameters);
  // std::cerr << " mesh->getDimension() " << mesh->getDimension() << "\n";
  // assert(0);
  assert(mesh);
  _alpha = parameters.doubles["alpha"];
  _diff = parameters.doubles["diff"];

  if(parameters.strings["beta"]=="zero")
  {
    _beta = std::unique_ptr<solvers::InitialConditionInterface>(new BetaZero());
  }
  else if(parameters.strings["beta"]=="east")
  {
    _beta = std::unique_ptr<solvers::InitialConditionInterface>(new BetaEast());
  }
  else
  {
    _beta = std::unique_ptr<solvers::InitialConditionInterface>(new BetaConstant());
  }
}

/*--------------------------------------------------------------------------*/
double Model::diffusion(double x, double y, double z)const
{
  return _diff;
}
void Model::beta(arma::vec& beta, double x, double y, double z, double t)const
{
  (*_beta)(beta, x, y, z, t);
}
void Model::reaction(arma::subview_col<double> f, const arma::subview_col<double> u)const
{
  f[0] = _alpha*u[0];
}
void Model::reaction_d(arma::mat& df, const arma::subview_col<double> u)const
{
  df(0,0) = _alpha;
}

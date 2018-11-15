#ifndef ___Model_hpp
#define ___Model_hpp

#include  "Solvers/modelinterface.hpp"

/*--------------------------------------------------------------------------*/
class Model : public solvers::ModelInterface
{
protected:
  std::shared_ptr<solvers::InitialConditionInterface> _beta;

public:
  ~Model();
  Model();
  Model( const Model& model);
  Model& operator=( const Model& model);
  std::string getClassName() const;
  std::string getInfo() const;

  double _alpha, _diff, _robin;
  void initModel(const mesh::MeshUnitInterface* mesh, const solvers::Parameters& parameters);

  double diffusion(double x, double y, double z)const;
  void beta(arma::vec& beta, double x, double y, double z, double t=0)const;
  void reaction(arma::subview_col<double> f, const arma::subview_col<double> u)const;
  void reaction_d(arma::mat& df, const arma::subview_col<double> u)const;
  std::shared_ptr<solvers::InitialConditionInterface> getBeta() const {return _beta;}
};

/*--------------------------------------------------------------------------*/
#endif

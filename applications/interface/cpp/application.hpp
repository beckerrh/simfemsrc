#ifndef ___Application_hpp
#define ___Application_hpp

#include  "Solvers/applicationinterface.hpp"
#include  "model.hpp"

/*--------------------------------------------------------------------------*/
class CutExactSolutionQuadraticCircle : public solvers::FunctionInterface
{
protected:
  mutable arma::vec _p;
  double _r02, _k1, _k2;
  const solvers::FunctionInterface* _phi;
public:
  CutExactSolutionQuadraticCircle(double K1, double k2);
  std::string getName() const;
  void setPhi(const solvers::FunctionInterface* phi);
  void operator()(arma::vec& u, double x, double y, double z, double t = 0.) const;
  void x (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void y (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void z (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void t (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void xx(arma::vec& u, double x, double y, double z, double t = 0.) const;
  void yy(arma::vec& u, double x, double y, double z, double t = 0.) const;
  void zz(arma::vec& u, double x, double y, double z, double t = 0.) const;
};
/*--------------------------------------------------------------------------*/
class CutExactSolutionLinearStraight : public solvers::FunctionInterface
{
protected:
  double _xgamma, _k1, _k2, _p1, _p2;
public:
  CutExactSolutionLinearStraight(double xgamma, double K1, double k2);
  std::string getName() const;
  void operator()(arma::vec& u, double x, double y, double z, double t = 0.) const;
  void x (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void y (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void z (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void t (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void xx(arma::vec& u, double x, double y, double z, double t = 0.) const;
  void yy(arma::vec& u, double x, double y, double z, double t = 0.) const;
  void zz(arma::vec& u, double x, double y, double z, double t = 0.) const;
};
/*--------------------------------------------------------------------------*/
class CutExactSolutionQuadraticStraight : public solvers::FunctionInterface
{
protected:
  double _xgamma, _k1, _k2;
public:
  CutExactSolutionQuadraticStraight(double xgamma, double K1, double k2);
  std::string getName() const;
  void operator()(arma::vec& u, double x, double y, double z, double t = 0.) const;
  void x (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void y (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void z (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void t (arma::vec& u, double x, double y, double z, double t = 0.) const;
  void xx(arma::vec& u, double x, double y, double z, double t = 0.) const;
  void yy(arma::vec& u, double x, double y, double z, double t = 0.) const;
  void zz(arma::vec& u, double x, double y, double z, double t = 0.) const;
};

/*--------------------------------------------------------------------------*/
class Application : public solvers::ApplicationInterface
{
protected:
  std::shared_ptr<solvers::FunctionInterface> newExactSolution(std::string varname) const;
  std::unique_ptr<solvers::RightHandSideInterface> newRightHandSide(std::string varname) const;
  std::unique_ptr<solvers::DirichletInterface> newDirichlet(std::string varname) const;
  std::shared_ptr<solvers::FunctionInterface> newDataFunction(std::string varname) const;
  std::string _applicationname;
  const Model* _localmodel;

public:
  ~Application();
  Application(std::string applicationname);
  Application( const Application& application);
  Application& operator=( const Application& application);
  std::string getClassName() const;

  void initApplication(const mesh::MeshUnitInterface* mesh, const alat::StringVector& varnames, const alat::StringVector& varnamesdata, const alat::armaivec& ncomps, const alat::armaivec& ncompsdata);
  bool isStrongDirichlet(int color)const;
};

/*--------------------------------------------------------------------------*/
#endif

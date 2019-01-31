#ifndef ___Application_hpp
#define ___Application_hpp

#include  "Solvers/applicationinterface.hpp"
#include  "model.hpp"

/*--------------------------------------------------------------------------*/
class RDcosh : public solvers::FunctionInterface
{
protected:
  double _b, _c;
public:
  RDcosh(double diff, double alpha);
  std::string getName() const;
  void operator()(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void x (alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void y (alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void z (alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void t (alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void xx(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void yy(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void zz(alat::armavec& u, double x, double y, double z, double t = 0.) const;
};
class CDExpLayer : public solvers::FunctionInterface
{
protected:
  double _diff;
public:
  CDExpLayer(double diff);
  std::string getName() const;
  void operator()(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void x (alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void y (alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void z (alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void t (alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void xx(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void yy(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  void zz(alat::armavec& u, double x, double y, double z, double t = 0.) const;
};

/*--------------------------------------------------------------------------*/
class Application : public solvers::ApplicationInterface
{
protected:
  const Model* _localmodel;
  std::string _exactsolution;
  bool _strongdir;

    std::shared_ptr<solvers::FunctionInterface> newExactSolution(std::string varname) const;
    std::unique_ptr<solvers::RightHandSideInterface> newRightHandSide(std::string varname) const;
    std::unique_ptr<solvers::DirichletInterface> newDirichlet(std::string varname) const;
    std::shared_ptr<solvers::FunctionInterface> newDataFunction(std::string varname) const;

public:
  ~Application();
  Application(std::string exactsolution="none", bool strongdir=true);
  Application( const Application& application);
  Application& operator=( const Application& application);
  std::string getClassName() const;

  void initApplication(const mesh::MeshUnitInterface* mesh, const alat::StringVector& varnames, const alat::StringVector& varnamesdata, const alat::armaivec& ncomps, const alat::armaivec& ncompsdata);
  bool isStrongDirichlet(int color)const;
};

/*--------------------------------------------------------------------------*/
#endif

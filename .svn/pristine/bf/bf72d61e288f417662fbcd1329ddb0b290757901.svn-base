#ifndef __Solvers_FunctionInterface_hpp
#define __Solvers_FunctionInterface_hpp

#include  "Alat/interfacebase.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  /*--------------------------------------------------------------------------*/
  class FunctionInterface: public alat::InterfaceBase
  {
  protected:
      int _dim, _ncomp;
      double _eps, _epsinv;
      mutable arma::vec _ul, _ur;

  public:
    ~FunctionInterface();
    FunctionInterface();
    FunctionInterface( const FunctionInterface& function);
    FunctionInterface& operator=( const FunctionInterface& function);
    virtual std::string getClassName()const;
    void setNcompAndDim(int nocomp, int dim);

    virtual void operator()(arma::vec& u, double xv, double yv, double zv, double t=0.0)const=0;
    virtual void d (arma::vec& u, int i, double xv, double yv, double zv, double t = 0.) const;
    virtual void dd(arma::vec& u, int i, int j, double xv, double yv, double zv, double t = 0.) const;
    virtual void x (arma::vec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void y (arma::vec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void z (arma::vec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void xx(arma::vec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void yy(arma::vec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void zz(arma::vec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void xy(arma::vec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void yz(arma::vec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void xz(arma::vec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void t (arma::vec& u, double xv, double yv, double zv, double t = 0.) const;
  };

  /*--------------------------------------------------------------------------*/
  class InitialConditionInterface: public FunctionInterface
  {
  public:
    virtual std::string getClassName()const {return "InitialConditionInterface";}
  };
  /*--------------------------------------------------------------------------*/
  class RightHandSideInterface: public FunctionInterface
  {
  public:
    virtual std::string getClassName()const {return "RightHandSideInterface";}
  };
  class RightHandSideConstant : public RightHandSideInterface
  {
  private:
    double _d;
  public:
    RightHandSideConstant(double d=0.0);
    std::string getClassName()const;
    void operator()(arma::vec& f, double x, double y, double z, double t)const;
  };
  /*--------------------------------------------------------------------------*/
  class DirichletInterface: public FunctionInterface
  {
  public:
    virtual std::string getClassName()const {return "DirichletInterface";}
  };
  class DirichletConstant: public solvers::DirichletInterface
  {
  private:
    double _d;
  public:
    DirichletConstant(double d=0.0);
    std::string getClassName()const;
    void operator()(arma::vec& u, double x, double y, double z, double t)const;
  };
  class DirichletExactSolution : public solvers::DirichletInterface{
  private:
    const FunctionInterface& solution;
  public:
    DirichletExactSolution(const FunctionInterface& solution_in);
    std::string getClassName()const;
    void operator()(arma::vec& u, double x, double y, double z, double t)const;
  };
  /*--------------------------------------------------------------------------*/
  class NeumannInterface: public FunctionInterface
  {
  public:
    virtual std::string getClassName()const {return "NeumannInterface";}
    void operator()(arma::vec& u, double xv, double yv, double zv, double t=0.0)const{assert(0);}
    virtual void operator()(arma::vec& u, double nx, double ny, double nz, double x, double y, double z, double t=0.0)const=0;
  };
  class NeumannZero: public NeumannInterface
  {
  public:
    std::string getClassName()const;
    void operator()(arma::vec& u, double nx, double ny, double nz, double x, double y, double z, double t=0.0)const;
  };

}

/*--------------------------------------------------------------------------*/
#endif

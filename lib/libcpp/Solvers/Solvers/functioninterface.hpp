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
      mutable alat::armavec _ul, _ur;

  public:
    ~FunctionInterface();
    FunctionInterface();
    FunctionInterface( const FunctionInterface& function);
    FunctionInterface& operator=( const FunctionInterface& function);
    virtual std::string getClassName()const;
    void setNcompAndDim(int nocomp, int dim);

    virtual void operator()(alat::armavec& u, double xv, double yv, double zv, double t=0.0)const=0;
    virtual void d (alat::armavec& u, int i, double xv, double yv, double zv, double t = 0.) const;
    virtual void dd(alat::armavec& u, int i, int j, double xv, double yv, double zv, double t = 0.) const;
    virtual void x (alat::armavec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void y (alat::armavec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void z (alat::armavec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void xx(alat::armavec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void yy(alat::armavec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void zz(alat::armavec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void xy(alat::armavec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void yz(alat::armavec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void xz(alat::armavec& u, double xv, double yv, double zv, double t = 0.) const;
    virtual void t (alat::armavec& u, double xv, double yv, double zv, double t = 0.) const;
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
    void operator()(alat::armavec& f, double x, double y, double z, double t)const;
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
    void operator()(alat::armavec& u, double x, double y, double z, double t)const;
  };
  class DirichletExactSolution : public solvers::DirichletInterface{
  private:
    const FunctionInterface& solution;
  public:
    DirichletExactSolution(const FunctionInterface& solution_in);
    std::string getClassName()const;
    void operator()(alat::armavec& u, double x, double y, double z, double t)const;
  };
  /*--------------------------------------------------------------------------*/
  class NeumannInterface: public FunctionInterface
  {
  public:
    virtual std::string getClassName()const {return "NeumannInterface";}
    void operator()(alat::armavec& u, double xv, double yv, double zv, double t=0.0)const{assert(0);}
    virtual void operator()(alat::armavec& u, double nx, double ny, double nz, double x, double y, double z, double t=0.0)const=0;
  };
  class NeumannZero: public NeumannInterface
  {
  public:
    std::string getClassName()const;
    void operator()(alat::armavec& u, double nx, double ny, double nz, double x, double y, double z, double t=0.0)const;
  };

}

/*--------------------------------------------------------------------------*/
#endif

#ifndef __Solvers_ExactSolutions_h
#define __Solvers_ExactSolutions_h

#include  "Solvers/functioninterface.hpp"

/*---------------------------------------------------------*/
namespace solvers
{
  class ConstantSolution : public solvers::FunctionInterface
  {
protected:
    double _d;
public:
    ConstantSolution(double d = 1.2);
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
  class LinearSolution : public solvers::FunctionInterface
  {
private:
    double _a, _b, _c, _d;
public:
    LinearSolution(double a = 1., double b = 1.0, double c = 0.4, double d = 1.7);
    std::string getName() const;
    void operator()(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void x (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void y (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void z (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void xx(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void yy(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void zz(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  };

  class QuadraticSolution : public solvers::FunctionInterface
  {
private:
    double _a, _b, _c, _d, _e, _f, _g, _h, _i;
public:
    QuadraticSolution(double a = 1.2, double b = 0.0, double c = 1.1, double d = 1.1, double e = 1.1, double f = 1.1, double g = 1.1, double h = 1.1, double i = 1.1);
    std::string getName() const;
    void operator()(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void x (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void y (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void z (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void xx(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void yy(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void zz(alat::armavec& u, double vx, double vy, double vz, double t = 0.) const;
  };
  class CubicSolution : public solvers::FunctionInterface
  {
private:
    double _a, _b, _c, _d, _e, _f, _g, _h, _i, _j, _k, _l, _m, _n, _o, _p, _r, _s, _t, _u;
public:
    CubicSolution(double a = 1.2, double b = 0.0, double c = 1.1, double d = 1.1, double e = 1.1, double f = 1.1, double g = 1.7, double h = 1.8, double i = 1.9, double j = 2.0, double k = 0.0, double l = 0.0, double m = 1.1, double n = 1.1, double o = 1.1, double p = 1.1, double r = 1.7, double s = 1.8, double t = 1.9, double u = 2.0);
    std::string getName() const;
    void operator()(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void x (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void y (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void z (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void xx(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void yy(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void zz(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  };
  class CosinusSolution : public solvers::FunctionInterface
  {
private:
    double _a, _b, _c, _d;
public:
    CosinusSolution(double a = 2.0, double b = 1.0, double c = 1., double d = 2.);
    std::string getName() const;
    void operator()(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void x (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void y (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void z (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void xx(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void yy(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void zz(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  };
  class ExponentialSolution : public solvers::FunctionInterface
  {
private:
    double _x0, _y0, _z0, _eps;
    mutable double _r, _rx, _ry, _rz;
    void _init(double x, double y, double z) const;
public:
    ExponentialSolution(double x0 = 0.5, double y0 = 0.5, double z0 = 0.5, double eps = 0.01);
    std::string getName() const;
    void operator()(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void x (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void y (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void z (alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void xx(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void yy(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void zz(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  };
  class LDomainSolution : public solvers::FunctionInterface
  {
protected:
    double _theta(double x, double y, double z) const;

public:
    std::string getName() const;
    void operator()(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void x(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void y(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  };

  class SlitDomainSolution : public LDomainSolution
  {
public:
    std::string getName() const;
    void operator()(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void x(alat::armavec& u, double x, double y, double z, double t = 0.) const;
    void y(alat::armavec& u, double x, double y, double z, double t = 0.) const;
  };
}

/*---------------------------------------------------------*/

#endif

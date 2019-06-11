#include  "Solvers/exactsolutions.hpp"
#include  <cmath>

using namespace solvers;
using namespace std;

/*---------------------------------------------------------*/
ConstantSolution::ConstantSolution(double d) : _d(d) {}
string ConstantSolution::getClassName() const{return "ConstantSolution";}
void ConstantSolution::operator()(alat::armavec& u, double x, double y, double z, double t) const{u[0] = _d;}
void ConstantSolution::x (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void ConstantSolution::y (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void ConstantSolution::z (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void ConstantSolution::t (alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void ConstantSolution::xx(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void ConstantSolution::yy(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void ConstantSolution::zz(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.;}

/*---------------------------------------------------------*/
LinearSolution::LinearSolution(double a, double b, double c, double d) : _a(a), _b(b), _c(c), _d(d) {}
string LinearSolution::getClassName() const{return "LinearSolution";}
void LinearSolution::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _a+_b*x+_c*y+_d*z;
}
void LinearSolution::x (alat::armavec& u, double x, double y, double z, double t) const{u[0] = _b;}
void LinearSolution::y (alat::armavec& u, double x, double y, double z, double t) const{u[0] = _c;}
void LinearSolution::z (alat::armavec& u, double x, double y, double z, double t) const{u[0] = _d;}
void LinearSolution::xx(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void LinearSolution::yy(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.0;}
void LinearSolution::zz(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 0.;}

/*---------------------------------------------------------*/
QuadraticSolution::QuadraticSolution(double a, double b, double c, double d, double e, double f, double g, double h, double i) : _a(a), _b(b), _c(c), _d(d), _e(e), _f(f), _g(g), _h(h), _i(i) {}
string QuadraticSolution::getClassName() const{return "QuadraticSolution";}
void QuadraticSolution::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _a+_b*x+_c*y+_d*x*x+_e*y*y+_f*x*y + _g*x*z + _h*y*z + _i*z*z;
}
void QuadraticSolution::x(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _b+2.0*_d*x+_f*y + _g*z;
}
void QuadraticSolution::y(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _c+2.0*_e*y+_f*x + _h*z;
}
void QuadraticSolution::z(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _g*x + _h*y + 2.0*_i*z;
}
void QuadraticSolution::xx(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 2.0*_d;}
void QuadraticSolution::yy(alat::armavec& u, double x, double y, double z, double t) const{u[0] = 2.0*_e;}
void QuadraticSolution::zz(alat::armavec& u, double vx, double vy, double vz, double t) const{u[0] = 2.0*_i;}

/*---------------------------------------------------------*/
CubicSolution::CubicSolution(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j, double k, double l, double m, double n, double o, double p, double r, double s, double t, double u) : _a(a), _b(b), _c(c), _d(d), _e(e), _f(f), _g(g), _h(h), _i(i), _j(j), _k(k), _l(l), _m(m), _n(n), _o(o), _p(p), _r(r), _s(s), _t(t), _u(u) {}

string CubicSolution::getClassName() const
{
  return "CubicSolution";
}

void CubicSolution::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _a+_b*x+_c*y+_d*x*x+_e*y*y+_f*x*y+_g*x*x*x+_h*x*x*y+_i*x*y*y+_j*y*y*y+_k*z+_l*z*z+_m*x*z+_n*y*z+_o*z*z*z+_p*x*x*z+_r*y*y*z+_s*x*y*z+_t*x*z*z+_u*y*z*z;
}

void CubicSolution::x(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _b+2.0*_d*x+_f*y+3.0*_g*x*x+2.0*_h*x*y+_i*y*y+_m*z+2.0*_p*x*z+_s*y*z+_t*z*z;
}

void CubicSolution::y(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _c+2.0*_e*y+_f*x+_h*x*x+2.0*_i*x*y+3.0*_j*y*y+_n*z+2.0*_r*y*z+_s*x*z+_u*z*z;
}

void CubicSolution::z(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _k+2.0*_l*z+_m*x+_n*y+3.0*_o*z*z+_p*x*x+_r*y*y+_s*x*y+2.0*_t*x*z+2.0*_u*y*z;
}

void CubicSolution::xx(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = 2.0*_d+6.0*_g*x+2.0*_h*y+2.0*_p*z;
}

void CubicSolution::yy(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = 2.0*_e+2.0*_i*x+6.0*_j*y+2.0*_r*z;
}

void CubicSolution::zz(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = 2.0*_l+6.0*_o*z+2.0*_t*x+2.0*_u*y;
}

/*---------------------------------------------------------*/
CosinusSolution::CosinusSolution(double a, double b, double c, double d) : _a(a), _b(b), _c(c), _d(d) {}
string CosinusSolution::getClassName() const
{
  return "CosinusSolution";
}
void CosinusSolution::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = _a*cos(_b*M_PI*x)*cos(_c*M_PI*y)*cos(_d*M_PI*z);
}
void CosinusSolution::x(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = -_a*_b*M_PI*sin(_b*M_PI*x)*cos(_c*M_PI*y)*cos(_d*M_PI*z);
}
void CosinusSolution::y(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = -_a*_c*M_PI*cos(_b*M_PI*x)*sin(_c*M_PI*y)*cos(_d*M_PI*z);
}
void CosinusSolution::z(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = -_a*_d*M_PI*cos(_b*M_PI*x)*cos(_c*M_PI*y)*sin(_d*M_PI*z);
}
void CosinusSolution::xx(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = -_a*_b*M_PI*_b*M_PI*cos(_b*M_PI*x)*cos(_c*M_PI*y)*cos(_d*M_PI*z);
}
void CosinusSolution::yy(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = -_a*_c*M_PI*_c*M_PI*cos(_b*M_PI*x)*cos(_c*M_PI*y)*cos(_d*M_PI*z);
}
void CosinusSolution::zz(alat::armavec& u, double x, double y, double z, double t) const
{
  u[0] = -_a*_d*M_PI*_d*M_PI*cos(_b*M_PI*x)*cos(_c*M_PI*y)*cos(_d*M_PI*z);
}
/*---------------------------------------------------------*/
void ExponentialSolution::_init(double x, double y, double z) const
{
  double xm = x-_x0;
  double ym = y-_y0;
  double zm = z-_z0;
  _r = 0.5*( xm*xm + ym*ym + zm*zm )/_eps;
  _rx = xm;
  _ry = ym;
  _rz = zm;
}
ExponentialSolution::ExponentialSolution(double x0, double y0, double z0, double eps) : solvers::FunctionInterface(), _x0(x0), _y0(y0), _z0(z0), _eps(eps) {}
string ExponentialSolution::getClassName() const
{
  return "ExponentialSolution";
}
void ExponentialSolution::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  _init(x, y, z);
  u[0] = exp(-_r);
}
void ExponentialSolution::x(alat::armavec& u, double x, double y, double z, double t) const
{
  _init(x, y, z);
  u[0] = -_rx* exp(-_r)/_eps;
}
void ExponentialSolution::y(alat::armavec& u, double x, double y, double z, double t) const
{
  _init(x, y, z);
  u[0] = -_ry* exp(-_r)/_eps;
}
void ExponentialSolution::z(alat::armavec& u, double x, double y, double z, double t) const
{
  _init(x, y, z);
  u[0] = -_rz* exp(-_r)/_eps;
}
void ExponentialSolution::xx(alat::armavec& u, double x, double y, double z, double t) const
{
  _init(x, y, z);
  u[0] = ( _rx*_rx/_eps-1.0 )*exp(-_r)/_eps;
}
void ExponentialSolution::yy(alat::armavec& u, double x, double y, double z, double t) const
{
  _init(x, y, z);
  u[0] = ( _ry*_ry/_eps-1.0 )*exp(-_r)/_eps;
}
void ExponentialSolution::zz(alat::armavec& u, double x, double y, double z, double t) const
{
  _init(x, y, z);
  u[0] = ( _rz*_rz/_eps-1.0 )*exp(-_r)/_eps;
}

/*---------------------------------------------------------*/
double LDomainSolution::_theta(double x, double y, double z) const
{
  if(y >= 0.0)
  {
    return atan2(y, x);
  }
  else
  {
    return 2.0*M_PI-atan2(y, x);
  }
}

string LDomainSolution::getClassName() const
{
  return "LDomainSolution";
}

void LDomainSolution::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  double r2 = x*x+y*y;
  u[0] = pow(r2, 1.0/3.0)*sin(2.0*_theta(x, y, z)/3.0);
}

void LDomainSolution::x(alat::armavec& u, double x, double y, double z, double t) const
{
  double r2 = x*x+y*y;
  u[0] = ( 2.0/3.0 )*pow(r2, -2.0/3.0)*( sin(2.0*_theta(x, y, z)/3.0)*x-cos(2.0*_theta(x, y, z)/3.0)*y );
}

void LDomainSolution::y(alat::armavec& u, double x, double y, double z, double t) const
{
  double r2 = x*x+y*y;
  u[0] = ( 2.0/3.0 )*pow(r2, -2.0/3.0)*( sin(2.0*_theta(x, y, z)/3.0)*y+cos(2.0*_theta(x, y, z)/3.0)*x );
}

/*---------------------------------------------------------*/

string SlitDomainSolution::getClassName() const
{
  return "SlitDomainSolution";
}

void SlitDomainSolution::operator()(alat::armavec& u, double x, double y, double z, double t) const
{
  double r2 = x*x+y*y;
  u[0] = pow(r2, 1.0/4.0)*sin(_theta(x, y, z)/2.0);
}

void SlitDomainSolution::x(alat::armavec& u, double x, double y, double z, double t) const
{
  double r2 = x*x+y*y;
  u[0] = ( 1.0/2.0 )*pow(r2, -3.0/4.0)*( sin(_theta(x, y, z)/2.0)*x-cos(_theta(x, y, z)/2.0)*y );
}

void SlitDomainSolution::y(alat::armavec& u, double x, double y, double z, double t) const
{
  double r2 = x*x+y*y;
  u[0] = ( 1.0/2.0 )*pow(r2, -3.0/4.0)*( sin(_theta(x, y, z)/2.0)*y+cos(_theta(x, y, z)/2.0)*x );
}

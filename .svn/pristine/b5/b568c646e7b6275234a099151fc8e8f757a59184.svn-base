#include  "Alat/cubicspline.hpp"
#include  "Alat/vector.hpp"
#include  <cassert>
#include  <limits>

using namespace alat;

/*--------------------------------------------------------------------------*/
CubicSpline::~CubicSpline() {}
CubicSpline::CubicSpline(){}
CubicSpline::CubicSpline( const CubicSpline& cubicspline)
{
  assert(0);
}
CubicSpline& CubicSpline::operator=( const CubicSpline& cubicspline)
{
  assert(0);
  return *this;
}
std::string CubicSpline::getClassName() const
{
  return "CubicSpline";
}

std::ostream& alat::operator<<(std::ostream& os, const CubicSplineNewtonInfo& info)
{
  os << info.g << " " << info.h << " " << info.d;
  return os;
}

/*--------------------------------------------------------------------------*/
void CubicSpline::bezier(arma::vec3& r, int i, double t) const
{
  double s = 1-t;
  r = s*s*s*_C(i,0) + 3*s*s*t*_C(i,1) + 3*s*t*t*_C(i,2) + t*t*t*_C(i,3);
}
void CubicSpline::dbezier(arma::vec3& dr, int i, double t) const
{
  double s = 1-t;
  dr = 3*s*s*(_C(i,1)-_C(i,0)) + 6*s*t*(_C(i,2)-_C(i,1)) +  3*t*t*(_C(i,3)-_C(i,2));
}
void CubicSpline::ddbezier(arma::vec3& ddr, int i, double t) const
{
  double s = 1-t;
  ddr = 6*s*(_C(i,0)-2*_C(i,1)+_C(i,2)) + 6*t*(_C(i,1)-2*_C(i,2)+_C(i,3));
}

/*--------------------------------------------------------------------------*/
void CubicSpline::_testsmooth() const
{
  int n = _C.n_rows;
  Node r1, r2;
  for(int i=0;i<n-1;i++)
  {
    dbezier(r1, i, 1.0);
    dbezier(r2, i+1, 0.0);
    assert(arma::norm(r1-r1)<1e-12);
    ddbezier(r1, i, 1.0);
    ddbezier(r2, i+1, 0.0);
    assert(arma::norm(r1-r1)<1e-12);
  }
  dbezier(r1, n-1, 1.0);
  dbezier(r2, 0, 0.0);
  assert(arma::norm(r1-r1)<1e-12);
  ddbezier(r1, n-1, 1.0);
  ddbezier(r2, 0, 0.0);
  assert(arma::norm(r1-r1)<1e-12);
}
/*--------------------------------------------------------------------------*/
void CubicSpline::_computeNewtonInfo(CubicSplineNewtonInfo& info, const arma::vec3& p, int i, double t) const
{
  bezier(info.r, i, t);
  info.r -= p;
  dbezier(info.dr, i, t);
  ddbezier(info.ddr, i, t);
  info.g = arma::dot(info.r,info.dr);
  info.h = arma::dot(info.dr,info.dr)+arma::dot(info.r,info.ddr);
  info.d = arma::norm(info.r);
}

/*--------------------------------------------------------------------------*/
double CubicSpline::_newtonInterval(double t0, double t1, CubicSplineNewtonInfo& info, const arma::vec3& p, int i) const
{
  double tmin=t0, tmax=t1;

  CubicSplineNewtonInfo infoL, infoM, infoR;
  _computeNewtonInfo(infoL, p, i, tmin);
  bool condL = infoL.g > 0.0;
  _computeNewtonInfo(infoR, p, i, tmax);
  bool condR = infoR.g < 0.0;
  double t = 0.5*(tmin+tmax);
  _computeNewtonInfo(infoM, p, i, t);
  // std::cerr << "condL="<<condL << " condR="<<condR<<" infoL.d="<<infoL.d<<" infoM.d="<<infoM.d<<" infoR.d="<<infoR.d <<"\n";

  if(infoL.d <= fmin(infoR.d, infoM.d))
  {
    info = infoL;
    tmax = t;
    t = 0.5*(tmin+tmax);
    if(condL) {return tmin;}
  }
  else if(infoR.d <= fmin(infoL.d, infoM.d))
  {
    info = infoR;
    tmin = t;
    t = 0.5*(tmin+tmax);
    if(condR) {return tmax;}
  }
  // std::cerr << "init t="<< t << "\n";
  int niter = 50;
  for(int iter=0;iter<niter;iter++)
  {
    _computeNewtonInfo(info, p, i, t);
    // std::cerr << iter << " NEWTON: " << t<< " ? " << info << "\n";
    if(abs(info.g)<1e-12) {return t;}
    if(info.h>0.0)
    {
      t -= info.g/info.h;
      if(t>tmax)
      {
        t = tmax;
        _computeNewtonInfo(info, p, i, t);
        if(info.g<0.0)
        {
          // std::cerr << iter << " QUITNEWTON: " << t<< " ? " << info << "\n";
          return t;
        }
      }
      else if(t<tmin)
      {
        t = tmin;
        _computeNewtonInfo(info, p, i, t);
        if(info.g>0.0)
        {
          // std::cerr << iter << " QUITNEWTON: " << t << " ? "<< info << "\n";
          return t;
        }
      }
    }
    else
    {
      double tN = t-info.g/arma::dot(info.dr,info.dr);
      tN = fmax(tmin,fmin(tN,tmax));
      _computeNewtonInfo(infoM, p, i, tN);
      std::cerr << "   t="<<t << " tN="<< tN << "\n";
      if(infoM.d<=info.d)
      {
        t = tN;
      }
      else if(info.g>0.0)
      {
        t = 0.5*(t+tmin);
      }
      else
      {
        t = 0.5*(t+tmax);
      }
    }
  }
  assert(0);
}

/*--------------------------------------------------------------------------*/
double CubicSpline::_newtonMinDist(CubicSplineNewtonInfo& info, const arma::vec3& p, int i) const
{
  int nint = 5;
  double dt = 1.0/(nint);
  alat::Vector<CubicSplineNewtonInfo> infos(5);
  alat::armavec ts(5);
  for(int iint=0;iint<nint;iint++)
  {
    ts[iint] = _newtonInterval(iint*dt, (iint+1)*dt, infos[iint], p, i);
  }
  double dmin=std::numeric_limits<float>::max();
  int imin = 0;
  for(int iint=0;iint<nint;iint++)
  {
    if(infos[iint].d < dmin)
    {
      dmin = infos[iint].d;
      imin =  iint;
    }
  }
  info = infos[imin];
  // std::cerr << "**** " << ts[imin] << " " << infos[imin] << "\n";
  return ts[imin];

  return _newtonInterval(0.0, 1.0, info, p, i);
  double t = 0.5;
  int niter = 50;
  for(int iter=0;iter<niter;iter++)
  {
    _computeNewtonInfo(info, p, i, t);
    std::cerr << iter << " NEWTON: " << t<< " ? " << info << "\n";
    if(abs(info.g)<1e-12) {return t;}
    if(info.h>0.0)
    {
      t -= info.g/info.h;
      if(t>1.0)
      {
        t = 1.0;
        _computeNewtonInfo(info, p, i, t);
        if(info.g<0.0)
        {
          std::cerr << iter << " QUITNEWTON: " << t<< " ? " << info << "\n";
          return t;
        }
      }
      else if(t<0.0)
      {
        t = 0.0;
        _computeNewtonInfo(info, p, i, t);
        if(info.g>0.0)
        {
          std::cerr << iter << " QUITNEWTON: " << t << " ? "<< info << "\n";
          return t;
        }
      }
    }
    else
    {
      if(info.g>0.0)
      {
        t = 0.5*t;
      }
      else
      {
        t = 1.0-2.0*t;
      }
    }
  }
  assert(0);
}

/*--------------------------------------------------------------------------*/
double CubicSpline::signeddistance(const arma::vec3& p) const
{
  int n = _C.n_rows;
  alat::Vector<CubicSplineNewtonInfo> infos(n);
  alat::armavec ts(n);
  for(int i=0;i<n;i++)
  {
    ts[i] =  _newtonMinDist(infos[i], p, i);
  }
  double dmin=std::numeric_limits<float>::max();
  int imin = 0;
  for(int i=0;i<n;i++)
  {
    if(infos[i].d < dmin)
    {
      dmin = infos[i].d;
      imin =  i;
    }
  }
  // std::cerr << "#### " << ts[imin] << " " << infos[imin] << "\n";
  if(infos[imin].r[1]*infos[imin].dr[0]-infos[imin].r[0]*infos[imin].dr[1] < 0)
  {
    return -1.0*dmin;
  }
  else
  {
    return dmin;
  }







  armavec dist(n);
  for(int i=0;i<n;i++)
  {
    dist[i] = arma::norm(_C(i,0)-p);
  }
  int i1 = dist.index_min();
  int i2 = i1-1;
  if(i1==0) {i2=n-1;}

  CubicSplineNewtonInfo info1, info2;
  double t1 = _newtonMinDist(info1, p, i1);
  double t2 = _newtonMinDist(info2, p, i2);

  double sign = 1.0;
  if(info1.d<=info2.d)
  {
    if(info1.r[1]*info1.dr[0]-info1.r[0]*info1.dr[1] < 0)
    {
      sign = -1.0;
    }
    return sign*info1.d;
  }
  else
  {
    if(info2.r[1]*info2.dr[0]-info2.r[0]*info2.dr[1]< 0)
    {
      sign = -1.0;
    }
    return sign*info2.d;
  }
}

/*--------------------------------------------------------------------------*/
void CubicSpline::generate(const armamat& points)
{
  // std::cerr << " points = " << points;
  assert(points.n_rows ==3);
  int n = points.n_cols;
  _C.set_size(n, 4);
  _A.resize(n,n);
  for(int i=0;i<n;i++)
  {
    _A(i,i) = 4.0;
  }
  for(int i=0;i<n-1;i++)
  {
    _A(i,i+1) = 1.0;
    _A(i+1,i) = 1.0;
  }
  _A(0,n-1) = 1.0;
  _A(n-1,0) = 1.0;
  _A = arma::inv(_A);
  arma::mat b(3,n, arma::fill::zeros);
  b.col(0) = 6*( points.col(1)-2*points.col(0)+points.col(n-1) );
  for(int i=1;i<n-1;i++)
  {
    b.col(i) = 6*(points.col(i-1)-2*points.col(i)+points.col(i+1));
  }
  b.col(n-1) = 6*(points.col(n-2)-2*points.col(n-1)+points.col(0));
  // std::cerr << " b = " << b;
  arma::mat M(3,n);
  M = b*_A;
  // std::cerr << " M = " << M;
  for(int i=0;i<n-1;i++)
  {
    _C(i, 0) = points.col(i);
    _C(i, 3) = points.col(i+1);
  }
  _C(n-1, 0) = points.col(n-1);
  _C(n-1, 3) = points.col(0);

  double unsurtrois = 1.0/3.0, unsurdixhuit = 1.0/18.0;
  for(int i=0;i<n-1;i++)
  {
    _C(i, 2) = unsurtrois*(points.col(i)+2*points.col(i+1)) - unsurdixhuit*(M.col(i)+2*M.col(i+1));
    _C(i, 1) = unsurtrois*(2*points.col(i)+points.col(i+1)) - unsurdixhuit*(2*M.col(i)+M.col(i+1));
  }
  _C(n-1, 2) = unsurtrois*(points.col(n-1)+2*points.col(0)) - unsurdixhuit*(M.col(n-1)+2*M.col(0));
  _C(n-1, 1) = unsurtrois*(2*points.col(n-1)+points.col(0)) - unsurdixhuit*(2*M.col(n-1)+M.col(0));
  _testsmooth();
  // for(int i=0;i<n;i++)
  // {
  //   std::cerr << i << "\n" << _C(i, 0).t() << _C(i, 1).t() << _C(i, 2).t() << _C(i, 3).t();
  // }
  // std::cerr << " _C = " << _C;
}

/*--------------------------------------------------------------------------*/
void CubicSpline::writeVtk(std::string filename) const
{
  int nc = _C.n_rows;
  int nn = 8*nc;
  alat::armamat newpoints(3,nn+1, arma::fill::zeros);

  int ic=0;
  double dtnc = 1.0/(double) (nc);
  double dtnn = 1.0/(double) (nn);
  arma::vec3 r;
  for(int i=0; i<nn;i++)
  {
    double t0 = i*dtnn;
    if(t0 > (ic+1)*dtnc) {ic++;}
    double t = nc*(t0-ic*dtnc);
    assert(ic < nc);
    bezier(r, ic, t);
    newpoints.col(i) = r;
  }
  newpoints.col(nn) = newpoints.col(0);

  std::ofstream file( filename.c_str());
  file<<"# vtk DataFile Version 4.0 "<<std::endl<<"output from SimFem"<<std::endl;
  file<<"ASCII"<<std::endl<<"DATASET POLYDATA"<<std::endl<<std::endl;

  file<<"POINTS "<<nn<<" FLOAT"<<std::endl;
  for(int i = 0; i < nn; i++)
  {
    file<<newpoints(0,i)<<" "<<newpoints(1,i)<<" "<<newpoints(2,i)<<" "<<std::endl;
  }
  file<<"\nLINES "<<nn<<" "<< 3*nn << std::endl;
  for(int i = 0; i < nn-1; i++)
  {
    file << 2<<" "<<i<<" "<<i+1<<" "<<std::endl;
  }
  file << 2<<" "<<nn-1<<" "<<0<<" "<<std::endl;
  file<<std::endl;

  file.close();
}

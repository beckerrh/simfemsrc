#include  "Alat/node.hpp"
#include  "FadalightMesh/curvedboundarydescriptionbase.hpp"
#include  <math.h>
#include  <iostream>
#include  <cassert>

using namespace FadalightMesh;
using namespace std;

/*---------------------------------------------------*/
CurvedBoundaryDescriptionBase::~CurvedBoundaryDescriptionBase()  {}
CurvedBoundaryDescriptionBase::CurvedBoundaryDescriptionBase() : FadalightMesh::CurvedBoundaryDescriptionInterface(), _parameters() {}
CurvedBoundaryDescriptionBase::CurvedBoundaryDescriptionBase(const CurvedBoundaryDescriptionBase& curvedboundarydescriptionbase) : FadalightMesh::CurvedBoundaryDescriptionInterface(curvedboundarydescriptionbase)
{
  _parameters.set_size(curvedboundarydescriptionbase._parameters.size());
  _parameters = curvedboundarydescriptionbase._parameters;
}

/*---------------------------------------------------*/

CurvedBoundaryDescriptionBase& CurvedBoundaryDescriptionBase::operator=(const CurvedBoundaryDescriptionBase& curvedboundarydescriptionbase)
{
  assert(0);
}

/*---------------------------------------------------*/

double CurvedBoundaryDescriptionBase::getParameter(int i) const
{
  return _parameters[i];
}

/*---------------------------------------------------*/

alat::armavec& CurvedBoundaryDescriptionBase::getParameters()
{
  return _parameters;
}

/*---------------------------------------------------*/

const alat::armavec& CurvedBoundaryDescriptionBase::getParameters() const
{
  return _parameters;
}

/*---------------------------------------------------*/

void CurvedBoundaryDescriptionBase::write(std::ostream& out, arma::file_type datatype) const
{
  _parameters.save(out, datatype);
}

/*---------------------------------------------------*/

void CurvedBoundaryDescriptionBase::read(std::istream& in)
{
  _parameters.load(in);
}

/*---------------------------------------------------*/

void CurvedBoundaryDescriptionBase::grad(alat::Node& dst, const alat::Node& src) const
{
  double eps = 1e-6;

  alat::Node cl(src), cr(src);
  for(int i = 0; i < 3; i++)
  {
    cl[i] -= eps;
    cr[i] += eps;
    dst[i] = ( ( *this )( cr.x(), cr.y(), cr.z() )-( *this )( cl.x(), cl.y(), cl.z() ) )/( 2.*eps );
    cl[i] += eps;
    cr[i] -= eps;
  }
}

/*---------------------------------------------------*/

double CurvedBoundaryDescriptionBase::newton(alat::Node& dst) const
{
  int maxi = 20;
  double rtol  = 1.e-14;

  alat::Node z;

  // grad(z, dst);
  // double res = ( *this )( dst.x(), dst.y(), dst.z() );
  double res = 2.*rtol, res0=1.0;

  for(int i = 0; ( i < maxi ) && ( fabs(res) > rtol*res0 ); i++)
  {
    alat::Node zz;
    grad(zz, dst);
    if(i==0)
    {
      z = zz;
      // std::cerr << "CurvedBoundaryDescriptionBase::newton() normal " << z << "\n";
    }
    double bgrad = z*zz;

    if(fabs(bgrad) <= 1.e-15)
    {
      cerr << "BoundaryDescrition::newton()\n";
      cerr << "Grad=0 in boundary_newton (res= "<< res << " )\n";
      cerr << "iter " << i << endl;
      cerr << "grad " << zz << endl;
      assert(0);
    }
    res /= bgrad;
    dst.add(-res, z);
    res = ( *this )( dst.x(), dst.y(), dst.z() );
    if(i==0)
    {
      res0 = fabs(res);
    }
    // cerr << "@@@ CurvedBoundaryDescriptionBase::newton() iter " << i << " " << res << endl;
  }
  if(fabs(res) > rtol)
  {
    cerr << "BoundaryDescrition::newton()\n";
    cerr << "No Convergence in boundary_newton (res= "<< res << " )\n";
    assert(0);
  }
  return fabs(res);
}

/*---------------------------------------------------*/

double CurvedBoundaryDescriptionBase::newton(alat::Node& dst, const alat::Node& x1, const alat::Node& x2) const
{
  int maxi = 20;
  double rtol  = 1.e-14;

  alat::Node z;

  double res = 2.*rtol, res0=1.0;

  for(int i = 0; ( i < maxi ) && ( fabs(res) > rtol*res0 ); i++)
  {
    alat::Node zz;
    grad(zz, dst);
    if(i==0)
    {
      z = zz;
      // std::cerr << "CurvedBoundaryDescriptionBase::newton() normal " << z << "\n";
    }
    double bgrad = z*zz;

    if(fabs(bgrad) <= 1.e-15)
    {
      cerr << "BoundaryDescrition::newton()\n";
      cerr << "Grad=0 in boundary_newton (res= "<< res << " )\n";
      cerr << "iter " << i << endl;
      assert(0);
    }
    res /= bgrad;
    dst.add(-res, z);
    res = ( *this )( dst.x(), dst.y(), dst.z() );
    if(i==0)
    {
      res0 = fabs(res);
    }
    // cerr << "@@@ CurvedBoundaryDescriptionBase::newton() iter " << i << " " << res << endl;
  }
  if(fabs(res) > rtol)
  {
    cerr << "BoundaryDescrition::newton()\n";
    cerr << "No Convergence in boundary_newton (res= "<< res << " )\n";
    assert(0);
  }
  return fabs(res);
}

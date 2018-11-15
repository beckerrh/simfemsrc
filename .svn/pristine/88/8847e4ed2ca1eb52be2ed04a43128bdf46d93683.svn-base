#ifndef  __Alat_VecOp_h
#define  __Alat_VecOp_h

#include  "Alat/armadillo.hpp"
#include  <cassert>

/*----------------------------------------------------------*/
namespace alat
{
  inline void vectorProduct(arma::vec& w, const arma::vec& u, const arma::vec&v)
  {
  	assert(w.size()==3);
  	assert(u.size()==3);
  	assert(v.size()==3);
    w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
  }
  // inline void vectorProduct(arma::subview_col<double> w, const arma::vec& u, const arma::vec&v)
  // {
  // 	assert(u.size()==3);
  // 	assert(v.size()==3);
  //   w[0] = u[1]*v[2] - u[2]*v[1];
  //   w[1] = u[2]*v[0] - u[0]*v[2];
  //   w[2] = u[0]*v[1] - u[1]*v[0];
  // }
}

#endif

#include  "Perulangan/fgmres.hpp"
#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  <cassert>
#include  <sstream>
#include  "Alat/armadillo.hpp"

using namespace perulangan;

/*--------------------------------------------------------------------------*/
FGmres::~FGmres(){}
FGmres::FGmres(int niteration) : perulangan::IterativeSolverWithPreconditioner(), _niteration(niteration){}
FGmres::FGmres( const FGmres& gmres) : IterativeSolverWithPreconditioner(gmres)
{
  assert(0);
}
FGmres& FGmres::operator=( const FGmres& gmres)
{
  IterativeSolverWithPreconditioner::operator=(gmres);
  assert(0);
  return *this;
}

std::string FGmres::getClassName() const
{
  std::stringstream ss;
  ss << _niteration;
  return "FGmres_"+ ss.str();
}

/*--------------------------------------------------------------------------*/

int FGmres::getNVectors() const
{
  return 2*( _niteration+1 );
}

/*--------------------------------------------------------------------------*/

void FGmres::solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const
{
  arma::mat _H;
  const perulangan::IterationInfo& info = *getIterationInfo();
  alat::GhostVector& r = getMemory(0);

  getVisitor()->residual(A, r, u, f);
  double betafirst = getVisitor()->vectorNorm(r);

  if( betafirst < info.getGlobalTol() )
  {
    status = perulanganEnums::IterationStatusNone;
    return;
  }
  info.checkIteration(status, betafirst, 0);

  getVisitor()->vectorScale(r, 1.0/betafirst);

  for(int iteration = 0; iteration < _niteration; iteration++)
  {
    int indexw = 2*( iteration+1 );
    getVisitor()->vectorZero( getMemory(2*iteration+1) );
    getPreconditioner()->solveApproximate( _statuspreconditioner, A, getMemory(2*iteration+1), getMemory(2*iteration), iteration );
    getVisitor()->vectorZero( getMemory(indexw) );
    getVisitor()->matrixVectorProduct(A, getMemory(indexw), getMemory(2*iteration+1), 1.0);
    _H.resize(iteration+2, iteration+1);
    for(int irow = 0; irow < iteration+1; irow++)
    {
      double scp = getVisitor()->vectorDot( getMemory(indexw), getMemory(2*irow) );
      getVisitor()->vectorAdd( getMemory(indexw), -scp, getMemory(2*irow) );
      _H(irow, iteration) = scp;
    }
    double beta = getVisitor()->vectorNorm( getMemory(indexw) );
    getVisitor()->vectorScale(getMemory(indexw), 1.0/beta);
    _H(iteration+1, iteration) = beta;
    // std::cerr << " _H = " << _H << "\n";
  }
  // arma::mat Q, R;
  // qr(Q, R, _H);
  // // qr_econ(Q,R,_H);
  // std::cerr << " Q = " << Q << "\n";
  // std::cerr << " R = " << R << "\n";

  arma::vec b(_niteration+1, arma::fill::zeros);
  b[0] = betafirst;
  // std::cerr << " b = " << b << "\n";
  arma::vec solution = arma::solve(_H, b);
  double res = arma::norm(b-_H*solution, 2);
  // std::cerr << " res = " << arma::norm(b-_H*solution, 2) << "\n";

  info.checkIteration(status, res, _niteration-1);
  bool reached = ((status==perulanganEnums::IterationStatusConverged) or (status==perulanganEnums::IterationStatusMaxIter));

  for(int i = 0; i < _niteration; i++)
  {
    getVisitor()->vectorAdd( u, solution[i], getMemory(2*i+1) );
  }
  if(reached)
  {
    status = perulanganEnums::IterationStatusConverged;
    return;
  }
  else
  {
    status = perulanganEnums::IterationStatusDiverged;
  }
}

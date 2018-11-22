#include  "Perulangan/bicgstab.hpp"
#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
BiCgStab::~BiCgStab(){}
BiCgStab::BiCgStab() : IterativeSolverWithPreconditioner(){}
BiCgStab::BiCgStab( const BiCgStab& bicgstab) : IterativeSolverWithPreconditioner(bicgstab)
{
  assert(0);
}

BiCgStab& BiCgStab::operator=( const BiCgStab& bicgstab)
{
  IterativeSolverWithPreconditioner::operator=(bicgstab);
  assert(0);
  return *this;
}

std::string BiCgStab::getClassName() const
{
  return "BiCgStab";
}

/*--------------------------------------------------------------------------*/
int BiCgStab::getNVectors() const
{
  return 7;
}

/*--------------------------------------------------------------------------*/


void BiCgStab::solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const
{
  const perulangan::IterationInfo& info = *getIterationInfo();
  double epsilon = 1e-32;

  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& rtilde = getMemory(1);
  alat::GhostVector& p = getMemory(2);
  alat::GhostVector& v = getMemory(3);
  alat::GhostVector& precvec = getMemory(4);
  alat::GhostVector& s = getMemory(5);
  alat::GhostVector& t = getMemory(6);

  getVisitor()->vectorEqual(r, f);
  getVisitor()->matrixVectorProduct(A, r, u, -1.0);
  double res = getVisitor()->vectorNorm(r);
  // std::cerr << "BiCgStab res " << res << "\n";
  if( res < info.getGlobalTol() )
  {
    status = perulanganEnums::IterationStatusNone;
    return;
  }
  getVisitor()->vectorEqual(rtilde, r);

  double rho, rhoold;
  double alpha, omega;

  bool reached = 0;
  for(int iteration = 0; !reached; iteration++)
  {
    rhoold = rho;
    rho = getVisitor()->vectorDot(r, rtilde);
    if(fabs(rho) < epsilon)
    {
      std::cerr << "¨*** ERROR in BiCgStab: rho is zero iteration rho " << iteration << " " << rho << "\n";
      assert(0);
    }
    if(iteration == 0)
    {
      getVisitor()->vectorEqual(p, r);
    }
    else
    {
      double beta = ( rho/rhoold )*( alpha/omega );
      getVisitor()->vectorAdd(p, beta-1.0, p);
      getVisitor()->vectorAdd(p, 1.0, r);
      getVisitor()->vectorAdd(p, -beta*omega, v);
    }
    // getVisitor()->precondition(precvec, p);
    getPreconditioner()->solveApproximate(_statuspreconditioner, A, precvec, p, iteration);
    getVisitor()->vectorZero(v);
    getVisitor()->matrixVectorProduct(A, v, precvec, 1.0);

    alpha = rho / getVisitor()->vectorDot(v, rtilde);
    getVisitor()->vectorEqual(s, r);
    getVisitor()->vectorAdd(s, -alpha, v);
    res = getVisitor()->vectorNorm(s);
    getVisitor()->vectorAdd(u, alpha, precvec);
    info.checkIteration(status, res, iteration);
    if( ( status == perulanganEnums::IterationStatusConverged ) or ( status == perulanganEnums::IterationStatusDiverged ) or (status == perulanganEnums::IterationStatusMaxIter) )
    {
      return;
    }
    // getVisitor()->precondition(precvec, s);
    getPreconditioner()->solveApproximate(_statuspreconditioner, A, precvec, s, iteration);
    getVisitor()->vectorZero(t);
    getVisitor()->matrixVectorProduct(A, t, precvec, 1.0);

    omega = getVisitor()->vectorDot(t, s)/getVisitor()->vectorDot(t, t);
    getVisitor()->vectorAdd(u, omega, precvec);
    getVisitor()->vectorEqual(r, s);
    getVisitor()->vectorAdd(r, -omega, t);
    if(fabs(omega) < epsilon)
    {
      std::cerr << "BiCgStab: omega is zero\n";
      assert(0);
    }
  }
  status = perulanganEnums::IterationStatusDiverged;
  return;
}

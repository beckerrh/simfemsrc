#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/orthominnew.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
OrthominNew::~OrthominNew(){}
OrthominNew::OrthominNew() : IterativeSolverWithPreconditioner(){}
OrthominNew::OrthominNew( const OrthominNew& orthomin) : IterativeSolverWithPreconditioner(orthomin)
{
  assert(0);
}
OrthominNew& OrthominNew::operator=( const OrthominNew& orthomin)
{
  IterativeSolverWithPreconditioner::operator=(orthomin);
  assert(0);
  return *this;
}
std::string OrthominNew::getClassName() const
{
  return "OrthominNew";
}
/*--------------------------------------------------------------------------*/

int OrthominNew::getNVectors() const
{
  return 7;
}

/*--------------------------------------------------------------------------*/

void OrthominNew::solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const
{
  const perulangan::IterationInfo& info = *getIterationInfo();

  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& z = getMemory(1);
  alat::GhostVector& p = getMemory(2);
  alat::GhostVector& y = getMemory(3);
  alat::GhostVector& Az = getMemory(4);
  alat::GhostVector& h = getMemory(5);
  alat::GhostVector& r2 = getMemory(6);

  getVisitor()->residual(A, r, u, f);
  getVisitor()->vectorEqual(r2, r);
  double res = getVisitor()->vectorNorm(r);
  if( res < info.getGlobalTol() )
  {
    status = perulanganEnums::IterationStatusNone;
    return;
  }
  int iteration=0;
  getPreconditioner()->solveApproximate(_statuspreconditioner, A, z, r, iteration);

  getVisitor()->vectorEqual(p, z);
  getVisitor()->vectorEqual(y, z);
  getVisitor()->vectorScale(y, -1.0);
  // double nu = getVisitor()->vectorDot(z,z);


  double omega, eta, nu;
  bool reached = 0;
  for(int iteration = 0; !reached; iteration++)
  {
    if( iteration == info.getMaxiter() )
    {
      status = perulanganEnums::IterationStatusDiverged;
      return;
    }
    getVisitor()->vectorZero(h);
    getVisitor()->matrixVectorProduct(A, h, z, 1.0);
    getPreconditioner()->solveApproximate(_statuspreconditioner, A, Az, h, iteration);

    double scpzAz = getVisitor()->vectorDot(Az,z);
    double scpAzAz = getVisitor()->vectorDot(Az,Az);
    if(iteration==0)
    {
      omega = scpzAz/scpAzAz;
      eta = 0.0;
      nu = omega*scpzAz;
    }
    else
    {
      double omegaold=omega;
      double scpyAz = getVisitor()->vectorDot(Az,y);
      double denominateur = nu*scpAzAz - scpyAz*scpyAz;
      omega = nu*scpzAz/denominateur;
      eta = -scpyAz*scpzAz/denominateur;
      getVisitor()->vectorScale(p, eta*omegaold/omega);
      getVisitor()->vectorAdd(p, 1.0, z);
    }
    getVisitor()->vectorAdd(u, omega, p);
    getVisitor()->vectorAdd(r2, eta, Az);
    getVisitor()->residual(A, r, u, f);
    res = getVisitor()->vectorNorm(r);
    // std::cerr << "res res2 " << res << " =? " << getVisitor()->vectorNorm(r2) << "\n";
    info.checkIteration(status, res, iteration);
    bool reached = status==perulanganEnums::IterationStatusConverged;
    if(reached)
    {
      return;
    }
    getVisitor()->vectorScale(y, eta);
    getVisitor()->vectorAdd(y, omega, Az);
    getVisitor()->vectorAdd(z, -1.0, y);

  }
  status = perulanganEnums::IterationStatusDiverged;
  return;
}

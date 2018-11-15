#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/orthomin.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
Orthomin::~Orthomin(){}
Orthomin::Orthomin() : IterativeSolverWithPreconditioner(){}
Orthomin::Orthomin( const Orthomin& orthomin) : IterativeSolverWithPreconditioner(orthomin)
{
  assert(0);
}
Orthomin& Orthomin::operator=( const Orthomin& orthomin)
{
  IterativeSolverWithPreconditioner::operator=(orthomin);
  assert(0);
  return *this;
}
std::string Orthomin::getClassName() const
{
  return "Orthomin";
}
/*--------------------------------------------------------------------------*/

int Orthomin::getNVectors() const
{
  return 5;
}

/*--------------------------------------------------------------------------*/

void Orthomin::solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const
{
  const perulangan::IterationInfo& info = *getIterationInfo();

  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& z = getMemory(1);
  alat::GhostVector& Az = getMemory(2);
  alat::GhostVector& zold = getMemory(3);
  alat::GhostVector& Azold = getMemory(4);

  getVisitor()->residual(A, r, u, f);
  double res = getVisitor()->vectorNorm(r);
  if( res < info.getGlobalTol() )
  {
    status = perulanganEnums::IterationStatusNone;
    return;
  }
  // getVisitor()->vectorZero(z);
  int iteration=0;
  getPreconditioner()->solveApproximate(_statuspreconditioner, A, z, r, iteration);
  bool reached = 0;
  for(int iteration = 0; !reached; iteration++)
  {
    if( iteration == info.getMaxiter() )
    {
      status = perulanganEnums::IterationStatusDiverged;
      return;
    }
    if(iteration)
    {
      getVisitor()->vectorEqual(Azold, Az);
    }
    getVisitor()->vectorZero(Az);
    getVisitor()->matrixVectorProduct(A, Az, z, 1.0);
    if(iteration)
    {
      double beta = getVisitor()->vectorDot(Az, Azold)/getVisitor()->vectorDot(Azold, Azold);
      getVisitor()->vectorAdd(z, -beta, zold);
      getVisitor()->vectorAdd(Az, -beta, Azold);
    }
    double scprAz = getVisitor()->vectorDot(r, Az);
    double normAz = getVisitor()->vectorDot(Az, Az);
    double omega = scprAz/normAz;
    // std::cerr << "iteration omega " << iteration << " " << omega << " scprAz normAz  " << scprAz << " " << normAz << "\n";
    getVisitor()->vectorAdd(u, omega, z);
    getVisitor()->vectorAdd(r, -omega, Az);
    res = getVisitor()->vectorNorm(r);
    // std::cerr << "Orthomin info res " << iteration << " " << res << "\n";
    info.checkIteration(status, res, iteration);
    bool reached = status==perulanganEnums::IterationStatusConverged;
    if(reached)
    {
      return;
    }
    getVisitor()->vectorEqual(zold, z);
    getVisitor()->vectorZero(z);
    getPreconditioner()->solveApproximate(_statuspreconditioner, A, z, r, iteration);
  }
  status = perulanganEnums::IterationStatusDiverged;
}

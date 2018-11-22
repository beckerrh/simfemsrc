#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  "Perulangan/richardsonmonotone.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/

RichardsonMonotone::~RichardsonMonotone()
{}

RichardsonMonotone::RichardsonMonotone() : IterativeSolverWithPreconditioner()
{}

RichardsonMonotone::RichardsonMonotone( const RichardsonMonotone& richardsonmonotone) : IterativeSolverWithPreconditioner(richardsonmonotone)
{
  assert(0);
}

RichardsonMonotone& RichardsonMonotone::operator=( const RichardsonMonotone& richardsonmonotone)
{
  IterativeSolverWithPreconditioner::operator=(richardsonmonotone);
  assert(0);
  return *this;
}

std::string RichardsonMonotone::getClassName() const
{
  return "RichardsonMonotone";
}

/*--------------------------------------------------------------------------*/

int RichardsonMonotone::getNVectors() const
{
  return 3;
}

/*--------------------------------------------------------------------------*/

void RichardsonMonotone::solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const
{
  const perulangan::IterationInfo& info = *getIterationInfo();

  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& z = getMemory(1);
  alat::GhostVector& Az = getMemory(2);

  bool reached = 0;
  double res;
  for(int iteration = 0; !reached; iteration++)
  {
    if( iteration == info.getMaxiter() )
    {
      status = perulanganEnums::IterationStatusDiverged;
      return;
    }
    if(iteration == 0)
    {
      getVisitor()->residual(A, r, u, f);
    }
    res = getVisitor()->vectorNorm(r);
    info.checkIteration(status, res, iteration);
    bool reached = status == perulanganEnums::IterationStatusConverged;
    if(reached)
    {
      return;
    }

    double rnorm = getVisitor()->vectorNorm(r);
    getVisitor()->vectorScale(r,1./rnorm);

    getVisitor()->vectorZero(z);
    getPreconditioner()->solveApproximate(_statuspreconditioner, A, z, r, iteration);



    getVisitor()->vectorZero(Az);
    getVisitor()->matrixVectorProduct(A, Az, z, 1.0);
    double scprAz = getVisitor()->vectorDot(r, Az);
    double normAz = getVisitor()->vectorDot(Az, Az);
    double omega = scprAz/normAz;
    // std::cerr << "iteration omega " << iteration << " " << omega << " scprAz normAz  " << scprAz << " " << normAz << " res=" << res << " omega=" << omega << "\n";
    if(normAz < 1e-28)
    {
      std::cerr << "iteration omega " << iteration << " " << omega << " scprAz normAz  " << scprAz << " " << normAz << " res=" << res << " omega=" << omega << "\n";
      omega = 1.0;
    }
    omega *= rnorm;
    getVisitor()->vectorAdd(u, omega, z);
    getVisitor()->vectorScale(r,rnorm);
    getVisitor()->vectorAdd(r, -omega, Az);
    // getVisitor()->residual(r, u, f);
  }
}

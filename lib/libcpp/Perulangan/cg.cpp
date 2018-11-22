#include  "Perulangan/cg.hpp"
#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
Cg::~Cg() {}
Cg::Cg(): IterativeSolverWithPreconditioner(){}
Cg::Cg( const Cg& cg): IterativeSolverWithPreconditioner(cg)
{
assert(0);
}

Cg& Cg::operator=( const Cg& cg)
{
IterativeSolverWithPreconditioner::operator=(cg);
assert(0);
return *this;
}

std::string Cg::getClassName() const
{
return "Cg";
}

/*--------------------------------------------------------------------------*/
int Cg::getNVectors() const
{
  return 4;
}

/*--------------------------------------------------------------------------*/
void Cg::solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const
{
  const perulangan::IterationInfo& info = *getIterationInfo();

  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& z = getMemory(1);
  alat::GhostVector& d = getMemory(2);
  alat::GhostVector& Ad = getMemory(3);

  getVisitor()->residual(A, r, u, f);
  double res = getVisitor()->vectorNorm(r);

  if( res < info.getGlobalTol() )
  {
    status = perulanganEnums::IterationStatusNone;
    return;
  }

  int iteration=0;
  getPreconditioner()->solveApproximate(_statuspreconditioner, A, z, r, iteration);
  double zr  = getVisitor()->vectorDot(z, r);
  getVisitor()->vectorEqual(d, z);

  bool reached = 0;

  for(int iteration = 0; !reached; iteration++)
  {
    getVisitor()->vectorZero(Ad);
    getVisitor()->matrixVectorProduct(A, Ad, d, 1.0);

    double alpha = zr / getVisitor()->vectorDot(d, Ad);

    getVisitor()->vectorAdd(r, -alpha, Ad);
    getVisitor()->vectorAdd(u, alpha, d);
    res = getVisitor()->vectorNorm(r);

    info.checkIteration(status, res, iteration);
    if( ( status == perulanganEnums::IterationStatusConverged ) or ( status == perulanganEnums::IterationStatusDiverged ) or (status == perulanganEnums::IterationStatusMaxIter) )
    {
      return;
    }

    getPreconditioner()->solveApproximate(_statuspreconditioner, A, z, r, iteration);

    double zrold = zr;
    zr   = getVisitor()->vectorDot(z, r);
    double beta = zr/zrold;

    getVisitor()->vectorAdd(d, beta-1.0, d);
    getVisitor()->vectorAdd(d, 1.0, z);
  }
  status = perulanganEnums::IterationStatusDiverged;
  return;
}

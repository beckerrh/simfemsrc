#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  "Perulangan/richardson.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/

Richardson::~Richardson()
{}

Richardson::Richardson(double omega) : IterativeSolverWithPreconditioner(), _omega(omega)
{}

Richardson::Richardson( const Richardson& richardson) : IterativeSolverWithPreconditioner(richardson)
{
  assert(0);
}

Richardson& Richardson::operator=( const Richardson& richardson)
{
  IterativeSolverWithPreconditioner::operator=(richardson);
  assert(0);
  return *this;
}

std::string Richardson::getClassName() const
{
  return "Richardson";
}

/*--------------------------------------------------------------------------*/

int Richardson::getNVectors() const
{
  return 2;
}

/*--------------------------------------------------------------------------*/

void Richardson::solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const
{
  const perulangan::IterationInfo& info = *getIterationInfo();

  alat::GhostVector& r = getMemory(0);
  alat::GhostVector& z = getMemory(1);

  bool reached = 0;
  for(int iteration = 0; !reached; iteration++)
  {
    if( iteration == info.getMaxiter() )
    {
      status = perulanganEnums::IterationStatusDiverged;
      return;
    }
    getVisitor()->residual(A, r, u, f);
    double res = getVisitor()->vectorNorm(r);
    info.checkIteration(status, res, iteration);
    bool reached = status==perulanganEnums::IterationStatusConverged;
    if(reached)
    {
      return;
    }
    // std::cerr << "_omega="<<_omega<<"\n";
    // getPreconditioner()->solveApproximate(status, A, u, f, iteration);
    getVisitor()->vectorZero(z);
    getPreconditioner()->solveApproximate(_statuspreconditioner, A, z, r, iteration);
    // getVisitor()->vectorZero(u);
    getVisitor()->vectorAdd(u, _omega, z);
  }
  status =  perulanganEnums::IterationStatusDiverged;
  return;
}

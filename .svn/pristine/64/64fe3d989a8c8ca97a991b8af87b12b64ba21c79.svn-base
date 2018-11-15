#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/trivialpreconditioner.hpp"
#include  <cassert>

using namespace perulangan;

/*--------------------------------------------------------------------------*/

TrivialPreconditioner::~TrivialPreconditioner()
{}

TrivialPreconditioner::TrivialPreconditioner() : perulangan::Preconditioner()
{}

TrivialPreconditioner::TrivialPreconditioner( const TrivialPreconditioner& trivialpreconditioner) : perulangan::Preconditioner(trivialpreconditioner)
{
  assert(0);
}

TrivialPreconditioner& TrivialPreconditioner::operator=( const TrivialPreconditioner& trivialpreconditioner)
{
  perulangan::Preconditioner::operator=(trivialpreconditioner);
  assert(0);
  return *this;
}

std::string TrivialPreconditioner::getClassName() const
{
  return "TrivialPreconditioner";
}

TrivialPreconditioner* TrivialPreconditioner::clone() const
{
  assert(0);
  return NULL;
//return new TrivialPreconditioner(*this);
}

/*--------------------------------------------------------------------------*/
void TrivialPreconditioner::setsmoothtype(std::string smoothtype)
{
  assert(0);
}

/*--------------------------------------------------------------------------*/

int TrivialPreconditioner::getNVectors() const
{
  return 0;
}

void TrivialPreconditioner::reInit()
{}
void TrivialPreconditioner::computePreconditioner()
{}
void TrivialPreconditioner::solveApproximate(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f, int iteration) const
{
  getVisitor()->vectorEqual(u, f);
}

#ifndef __Perulangan_RichardsonMonotone_h
#define __Perulangan_RichardsonMonotone_h

#include  "iterativesolverwithpreconditioner.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
class RichardsonMonotone : public IterativeSolverWithPreconditioner
{
private:
protected:
public:
  ~RichardsonMonotone();
  RichardsonMonotone();
  RichardsonMonotone( const RichardsonMonotone& richardsonmonotone);
  RichardsonMonotone& operator=( const RichardsonMonotone& richardsonmonotone);
  std::string getClassName() const;

  int getNVectors() const;
  void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
};
}

/*--------------------------------------------------------------------------*/

#endif

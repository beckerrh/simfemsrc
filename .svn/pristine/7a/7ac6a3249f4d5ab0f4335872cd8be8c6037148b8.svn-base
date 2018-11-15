#ifndef __Perulangan_BiCgStab_h
#define __Perulangan_BiCgStab_h

#include  "iterativesolverwithpreconditioner.hpp"

/*--------------------------------------------------------------------------*/
namespace perulangan
{
  class BiCgStab : public IterativeSolverWithPreconditioner
  {
public:
    ~BiCgStab();
    BiCgStab();
    BiCgStab( const BiCgStab& bicgstab);
    BiCgStab& operator=( const BiCgStab& bicgstab);
    std::string getClassName() const;

    int getNVectors() const;
    void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

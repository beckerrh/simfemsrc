#ifndef __Perulangan_OrthominNew_h
#define __Perulangan_OrthominNew_h

#include  "iterativesolverwithpreconditioner.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class OrthominNew : public IterativeSolverWithPreconditioner
  {
public:
    ~OrthominNew();
    OrthominNew();
    OrthominNew( const OrthominNew& orthomin);
    OrthominNew& operator=( const OrthominNew& orthomin);
    std::string getClassName() const;

    int getNVectors() const;
    void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

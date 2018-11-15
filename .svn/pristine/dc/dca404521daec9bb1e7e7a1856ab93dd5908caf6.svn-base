#ifndef __Perulangan_Orthomin_h
#define __Perulangan_Orthomin_h

#include  "iterativesolverwithpreconditioner.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class Orthomin : public IterativeSolverWithPreconditioner
  {
public:
    ~Orthomin();
    Orthomin();
    Orthomin( const Orthomin& orthomin);
    Orthomin& operator=( const Orthomin& orthomin);
    std::string getClassName() const;

    int getNVectors() const;
    void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

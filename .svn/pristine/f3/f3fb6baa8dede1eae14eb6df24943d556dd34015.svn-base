#ifndef __Perulangan_Cg_h
#define __Perulangan_Cg_h

#include  "iterativesolverwithpreconditioner.hpp"

/*--------------------------------------------------------------------------*/
namespace perulangan
{
  class Cg : public IterativeSolverWithPreconditioner
  {
  public:
    ~Cg();
    Cg();
    Cg( const Cg& cg);
    Cg& operator=( const Cg& cg);
    std::string getClassName() const;

    int getNVectors() const;
    void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

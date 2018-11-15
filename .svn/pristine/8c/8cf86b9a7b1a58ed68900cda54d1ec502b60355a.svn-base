#ifndef __Perulangan_FGmres_h
#define __Perulangan_FGmres_h

#include  "iterativesolverwithpreconditioner.hpp"

/*--------------------------------------------------------------------------*/
namespace perulangan
{
  class FGmres : public IterativeSolverWithPreconditioner
  {
protected:
    int _niteration;

public:
    ~FGmres();
    FGmres(int niteration);
    FGmres( const FGmres& gmres);
    FGmres& operator=( const FGmres& gmres);
    std::string getClassName() const;

    int getNVectors() const;
    void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

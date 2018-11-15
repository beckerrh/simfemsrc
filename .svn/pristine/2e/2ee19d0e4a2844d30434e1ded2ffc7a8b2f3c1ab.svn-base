#ifndef __Perulangan_Richardson_h
#define __Perulangan_Richardson_h

#include  "iterativesolverwithpreconditioner.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class Richardson : public IterativeSolverWithPreconditioner
  {
private:
protected:
    double _omega;
public:
    ~Richardson();
    Richardson(double omega = 1.0);
    Richardson( const Richardson& richardson);
    Richardson& operator=( const Richardson& richardson);
    std::string getClassName() const;

    int getNVectors() const;
    void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

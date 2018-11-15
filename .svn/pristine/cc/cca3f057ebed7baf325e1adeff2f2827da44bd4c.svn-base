#ifndef __Perulangan_SimpleIterativeSolver_h
#define __Perulangan_SimpleIterativeSolver_h

#include  "iterativesolverwithpreconditioner.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class SimpleIterativeSolver : public IterativeSolverWithPreconditioner
  {
private:
protected:
public:
    ~SimpleIterativeSolver();
    SimpleIterativeSolver();
    SimpleIterativeSolver( const SimpleIterativeSolver& simpleiterativesolver);
    SimpleIterativeSolver& operator=( const SimpleIterativeSolver& simpleiterativesolver);
    std::string getClassName() const;
    SimpleIterativeSolver* clone() const;

    void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

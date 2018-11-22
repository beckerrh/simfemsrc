#ifndef __Perulangan_NewtonSimple_h
#define __Perulangan_NewtonSimple_h

#include  "nonlinearsolver.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class NewtonSimple : public NonlinearSolver
  {
protected:
public:
    ~NewtonSimple();
    NewtonSimple();
    NewtonSimple( const NewtonSimple& newtonsimple);
    NewtonSimple& operator=( const NewtonSimple& newtonsimple);
    std::string getClassName() const;
    NewtonSimple* clone() const;

    int getNVectors() const;
    void solve(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f);
  };
}

/*--------------------------------------------------------------------------*/

#endif

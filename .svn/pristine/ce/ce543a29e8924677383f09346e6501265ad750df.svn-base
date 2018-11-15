#ifndef __Perulangan_NewtonRb_h
#define __Perulangan_NewtonRb_h

#include  "nonlinearsolver.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class NewtonRb : public NonlinearSolver
  {
protected:
    std::string _linesearchtype;

public:
    ~NewtonRb();
    NewtonRb(std::string linesearchtype);
    NewtonRb( const NewtonRb& newtonsimple);
    NewtonRb& operator=( const NewtonRb& newtonsimple);
    std::string getClassName() const;
    NewtonRb* clone() const;

    int getNVectors() const;
    void solve(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f);
  };
}

/*--------------------------------------------------------------------------*/

#endif

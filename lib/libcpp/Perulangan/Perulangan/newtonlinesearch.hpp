#ifndef __Perulangan_NewtonLineSearch_h
#define __Perulangan_NewtonLineSearch_h

#include  "nonlinearsolver.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class NewtonLineSearch : public NonlinearSolver
  {
protected:
    std::string _linesearchtype;

public:
    ~NewtonLineSearch();
    NewtonLineSearch(std::string linesearchtype);
    NewtonLineSearch( const NewtonLineSearch& newtonsimple);
    NewtonLineSearch& operator=( const NewtonLineSearch& newtonsimple);
    std::string getClassName() const;
    NewtonLineSearch* clone() const;

    int getNVectors() const;
    void solve(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f);
  };
}

/*--------------------------------------------------------------------------*/

#endif

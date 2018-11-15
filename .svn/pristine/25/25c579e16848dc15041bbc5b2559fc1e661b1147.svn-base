#ifndef __Perulangan_NewtonLavrentievStab_h
#define __Perulangan_NewtonLavrentievStab_h

#include  "newtonlavrentiev.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class NewtonLavrentievStab : public NewtonLavrentiev
  {
protected:
    bool _changeRegularization(int iteration, double& regularization, double modelfactor);

public:
    ~NewtonLavrentievStab();
    NewtonLavrentievStab();
    NewtonLavrentievStab( const NewtonLavrentievStab& newtonsimple);
    NewtonLavrentievStab& operator=( const NewtonLavrentievStab& newtonsimple);
    std::string getClassName() const;
    NewtonLavrentievStab* clone() const;

    int getNVectors() const;
    void solve(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f);
  };
}

/*--------------------------------------------------------------------------*/

#endif

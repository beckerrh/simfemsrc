#ifndef __Perulangan_NewtonLavrentiev_h
#define __Perulangan_NewtonLavrentiev_h

#include  "nonlinearsolver.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class NewtonLavrentiev : public NonlinearSolver
  {
protected:
    double _firstregularization, _maxregularization;
    std::string _type;

    bool _LINEARIZATION_TEST;
    bool _USE_LINEARIZATION;
    bool _changeRegularization(int iteration, double& regularization, double modelfactor);

public:
    ~NewtonLavrentiev();
    NewtonLavrentiev(std::string type="lin");
    NewtonLavrentiev( const NewtonLavrentiev& newtonsimple);
    NewtonLavrentiev& operator=( const NewtonLavrentiev& newtonsimple);
    std::string getClassName() const;
    NewtonLavrentiev* clone() const;

    // void basicInit(const alat::ParameterFile* parameterfile, std::string blockname);
    int getNVectors() const;
    void solve(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f);
  };
}

/*--------------------------------------------------------------------------*/

#endif

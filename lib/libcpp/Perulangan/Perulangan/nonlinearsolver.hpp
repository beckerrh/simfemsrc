#ifndef __Perulangan_NonlinearSolver_h
#define __Perulangan_NonlinearSolver_h

#include  "Alat/vector.hpp"
#include  "Alat/ghostvector.hpp"
#include  "iterationinfo.hpp"
#include  "nonlinearsolverinterface.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class NonlinearSolverVisitorInterface;

  struct NewtonInternalData
  {
    perulanganEnums::iterationstatus status_iteration;
    bool newmatrix, use_linesearch, use_regularization;
    int nliniter, iteration, irelax;
    double lintol, dxnorm, omega,residual, residual_old, modelfactor, regularization;
    NewtonInternalData();
  };
  std::ostream& operator<<(std::ostream& os, const NewtonInternalData& nd);

  class NonlinearSolver : public NonlinearSolverInterface
  {
private:
    int _id;
    static int _totalids;
    mutable alat::Vector<alat::GhostVector> _memory;
    std::shared_ptr<perulangan::NonlinearSolverVisitorInterface> _visitor;
    void _setNewtonInputData();

protected:
    NewtonInputData _newtoninputdata;
    mutable NewtonInternalData _newtoninternaldata;
    perulangan::IterationInfo _info;

    void memory();
    alat::GhostVector& getMemory(int i) const;
    bool _basicinitcalled;

    void _lineSearchMonotonicity(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostVector& u, const alat::GhostVector& du, alat::GhostVector& r, const alat::GhostVector& f) const;
    void _lineSearchArmijo(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostVector& u, const alat::GhostVector& du, alat::GhostVector& r, const alat::GhostVector& f, double kappa) const;
    void _lineSearchWolfe(perulangan::NewtonOutputData& newtonoutputdata, alat::GhostVector& u, const alat::GhostVector& du, alat::GhostVector& r, const alat::GhostVector& f, double kappa) const;
    bool _checkIteration(perulangan::NewtonOutputData& newtonoutputdata) const;

public:
    ~NonlinearSolver();
    NonlinearSolver();
    NonlinearSolver( const NonlinearSolver& nonlinearsolver);
    NonlinearSolver& operator=( const NonlinearSolver& nonlinearsolver);
    std::string getClassName() const;
    NonlinearSolver* clone() const;

    const NewtonInputData getNewtonInputData() const;
    void setNewtonInputData(const NewtonInputData newtoninputdata);
    void setVisitorPointer(std::shared_ptr<perulangan::NonlinearSolverVisitorInterface> visitor);
    const perulangan::NonlinearSolverVisitorInterface* getVisitor() const;
    perulangan::NonlinearSolverVisitorInterface* getVisitor();

    void init();
    std::ostream& printLoopInformation(std::ostream& os) const;
    void setVisitorPointer(std::unique_ptr<perulangan::NonlinearSolverVisitorInterface> visitor);
  };
}

/*--------------------------------------------------------------------------*/

#endif

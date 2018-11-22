#ifndef __Perulangan_NonlinearSolverInterface_h
#define __Perulangan_NonlinearSolverInterface_h

#include  "Alat/interfacebase.hpp"
#include  "enums.hpp"
#include  <memory>

/*--------------------------------------------------------------------------*/
namespace alat
{
  class GhostLinearSolver;
  class GhostMatrix;
  class GhostVector;
  class IntVector;
  class ParameterFile;
}
namespace perulangan
{
  class IterationInfo;
  class NewtonData;
  class NonlinearSolverVisitorInterface;
}
namespace perulangan
{
  struct NewtonInputData
  {
    bool printlinesearch;
    int maxnlinesearch, maxiter, printstep;
    double gtol, rtol, omegalinesearch, monotonyfactor, lineartoleranceincrease, rhomatrix;
    NewtonInputData();
    NewtonInputData(const NewtonInputData& newtoninputdata);
  };
  std::ostream& operator<<(std::ostream& os, const NewtonInputData& nd);
  struct NewtonOutputData
  {
    int niter_linear, nredo_matrix, niter;
    perulanganEnums::residualstatus residual_status;
    perulanganEnums::iterationstatus linear_solver_status;
    perulanganEnums::newtonstatus newton_status;
    perulanganEnums::matrixstatus matrix_status;
    NewtonOutputData();
  };
  std::ostream& operator<<(std::ostream& os, const NewtonOutputData& n);

  class NonlinearSolverInterface : public virtual alat::InterfaceBase
  {
  protected:
    std::string getInterfaceName() const;

  public:
    ~NonlinearSolverInterface();
    NonlinearSolverInterface();
    NonlinearSolverInterface( const NonlinearSolverInterface& linearsolverinterface);
    NonlinearSolverInterface& operator=( const NonlinearSolverInterface& linearsolverinterface);

    virtual void setVisitorPointer(std::shared_ptr<perulangan::NonlinearSolverVisitorInterface> visitor)=0;

    virtual const perulangan::NonlinearSolverVisitorInterface* getVisitor() const=0;
    virtual perulangan::NonlinearSolverVisitorInterface* getVisitor()=0;
    virtual void setNewtonInputData(const NewtonInputData newtoninputdata)=0;
    virtual const NewtonInputData getNewtonInputData() const=0;

    virtual int getNVectors() const;
    virtual void reInit();
    virtual void solve(perulangan::NewtonOutputData& newtondata, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f)=0;
    virtual std::ostream& printLoopInformation(std::ostream& os) const;
    virtual void init()=0;
    virtual std::ostream& write(std::ostream& os) const;
    virtual void addUpdate(perulanganEnums::iterationstatus& status, const alat::GhostVector& w, alat::GhostVector& u) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

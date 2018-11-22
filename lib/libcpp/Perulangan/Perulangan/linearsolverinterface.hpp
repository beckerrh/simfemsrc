#ifndef __Perulangan_LinearSolverInterface_h
#define __Perulangan_LinearSolverInterface_h

#include  "Perulangan/enums.hpp"
#include  "Alat/interfacebase.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class GhostMatrix;
  class GhostVector;
  class IntVector;
}

namespace perulangan
{
  class IterationInfo;

  class LinearSolverInterface : public virtual alat::InterfaceBase
  {
  protected:
    std::string getInterfaceName() const;

  public:
    ~LinearSolverInterface();
    LinearSolverInterface();
    LinearSolverInterface( const LinearSolverInterface& linearsolverinterface);
    LinearSolverInterface& operator=( const LinearSolverInterface& linearsolverinterface);

    virtual int getNVectors() const;
    virtual void reInit()=0;
    virtual void compute()=0;
    virtual void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const=0;
    virtual std::ostream& printLoopInformation(std::ostream& os) const;
    // virtual void basicInit(const alat::ParameterFile* parameterfile, std::string blockname);
    virtual std::ostream& write(std::ostream& os) const;
    virtual const IterationInfo* getIterationInfo() const;
    virtual IterationInfo* getIterationInfo();
    virtual void addUpdate(perulanganEnums::iterationstatus& status, const alat::GhostVector& w, alat::GhostVector& u) const;
    virtual bool hasIterationInfo() const = 0;
    virtual void restart();
    virtual void setTolerance(double rtol, double gtol);
    virtual int getNumberOfIterations() const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

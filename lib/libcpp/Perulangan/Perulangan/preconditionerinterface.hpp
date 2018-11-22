#ifndef __Perulangan_PreconditionerInterface_h
#define __Perulangan_PreconditionerInterface_h

#include  "enums.hpp"
#include  "Alat/interfacebase.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class GhostMatrix;
  class GhostVector;
  class IntVector;
  class ParameterFile;
  class StringVector;
}
namespace perulangan
{
  class IterativeSolverVisitorInterface;

  class PreconditionerInterface : public virtual alat::InterfaceBase
  {
protected:
    std::string getInterfaceName() const;

public:
    ~PreconditionerInterface();
    PreconditionerInterface();
    PreconditionerInterface( const PreconditionerInterface& vectorinterface);
    PreconditionerInterface& operator=( const PreconditionerInterface& vectorinterface);
    std::string getClassName() const;
    virtual std::ostream& printLoopInformation(std::ostream& os) const;

    virtual perulangan::IterativeSolverVisitorInterface* getVisitor();
    virtual const perulangan::IterativeSolverVisitorInterface* getVisitor() const;
    virtual int getNVectors() const = 0;
    // virtual void basicInit(const alat::ParameterFile* parameterfile, std::string blockname, perulangan::IterativeSolverVisitorInterface* visitor);
    virtual void reInit();
    virtual void computePreconditioner();
    virtual void solveApproximate(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f, int iteration) const;
    virtual std::ostream& write(std::ostream& os) const;
    virtual void fillzeros() const;
    virtual void setsmoothtype(std::string smoothtype) = 0;
  };
}

/*--------------------------------------------------------------------------*/

#endif

#ifndef __Perulangan_IterativeSolverWithPreconditioner_h
#define __Perulangan_IterativeSolverWithPreconditioner_h

#include  "iterativesolverwithvisitor.hpp"

/*--------------------------------------------------------------------------*/
namespace perulangan
{
  class PreconditionerInterface;

  class IterativeSolverWithPreconditioner : public virtual IterativeSolverWithVisitor
  {
private:
    PreconditionerInterface* _preconditioner;
    bool _deletepreconditioner;

protected:
    mutable perulanganEnums::iterationstatus _statuspreconditioner;

public:
    ~IterativeSolverWithPreconditioner();
    IterativeSolverWithPreconditioner();
    IterativeSolverWithPreconditioner( const IterativeSolverWithPreconditioner& iterativesolverwithpreconditioner);
    IterativeSolverWithPreconditioner& operator=( const IterativeSolverWithPreconditioner& iterativesolverwithpreconditioner);
    std::string getClassName() const;
    perulangan::IterativeSolverWithPreconditioner*  clone() const;

    // void basicInit(const alat::ParameterFile* parameterfile, std::string blockname);
    std::ostream& printLoopInformation(std::ostream& os) const;
    void reInit();
    void compute();

    perulangan::PreconditionerInterface*& newPreconditionerPointer();
    void setPreconditionerPointer(perulangan::PreconditionerInterface* preconditioner);
    const perulangan::PreconditionerInterface* getPreconditioner() const;
    perulangan::PreconditionerInterface* getPreconditioner();
    void addUpdate(perulanganEnums::iterationstatus& status, const alat::GhostVector& w, alat::GhostVector& u) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

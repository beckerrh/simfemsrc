#ifndef __Perulangan_Preconditioner_h
#define __Perulangan_Preconditioner_h

#include  "Alat/vector.hpp"
#include  "Alat/ghostvector.hpp"
#include  "preconditionerinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace perulangan
{
  class Preconditioner : public virtual PreconditionerInterface
  {
private:
    int _id;
    static int _totalids;

protected:
    mutable alat::Vector<alat::GhostVector> _memory;
    perulangan::IterativeSolverVisitorInterface* _visitor;
    std::string _smoothtype;
    alat::GhostVector& getMemory(int i) const;

public:
    ~Preconditioner();
    Preconditioner();
    Preconditioner( const Preconditioner& preconditioner);
    Preconditioner& operator=( const Preconditioner& preconditioner);
    std::string getClassName() const;
    Preconditioner* clone() const;

    perulangan::IterativeSolverVisitorInterface* getVisitor();
    const perulangan::IterativeSolverVisitorInterface* getVisitor() const;
    int getNVectors() const;
    // void basicInit(const alat::ParameterFile* parameterfile, std::string blockname, perulangan::IterativeSolverVisitorInterface* visitor);
    void memory();
    void setsmoothtype(std::string smoothtype);
  };
}

/*--------------------------------------------------------------------------*/

#endif

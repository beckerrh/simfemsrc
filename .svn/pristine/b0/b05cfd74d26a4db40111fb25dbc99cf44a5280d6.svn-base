#ifndef __Perulangan_IterativeSolverWithVisitor_h
#define __Perulangan_IterativeSolverWithVisitor_h

#include  "Alat/vector.hpp"
#include  "Alat/ghostvector.hpp"
#include  "iterationinfo.hpp"
#include  "linearsolverinterface.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class IterativeSolverVisitorInterface;

  class IterativeSolverWithVisitor : public virtual LinearSolverInterface
  {
private:
    int _id;
    static int _totalids;
    mutable alat::Vector<alat::GhostVector> _memory;
    perulangan::IterationInfo _info;
    perulangan::IterativeSolverVisitorInterface* _visitor;
    bool _deletevisitor;

protected:
    void memory();
    alat::GhostVector& getMemory(int i) const;
    bool _basicinitcalled;

public:
    ~IterativeSolverWithVisitor();
    IterativeSolverWithVisitor();
    IterativeSolverWithVisitor( const IterativeSolverWithVisitor& iterativesolverwithvisitor);
    IterativeSolverWithVisitor& operator=( const IterativeSolverWithVisitor& iterativesolverwithvisitor);
    std::string getClassName() const;
    perulangan::IterativeSolverWithVisitor* clone() const;

    // void basicInit(const alat::ParameterFile* parameterfile, std::string blockname);
    std::ostream& printLoopInformation(std::ostream& os) const;
    perulangan::IterativeSolverVisitorInterface*& newVisitorPointer();
    void setVisitorPointer(perulangan::IterativeSolverVisitorInterface* visitor);
    const perulangan::IterativeSolverVisitorInterface* getVisitor() const;
    perulangan::IterativeSolverVisitorInterface* getVisitor();
    const IterationInfo* getIterationInfo() const;
    IterationInfo* getIterationInfo();
    bool hasIterationInfo() const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

#ifndef __Perulangan_IterativeSolverVisitorInterface_h
#define __Perulangan_IterativeSolverVisitorInterface_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class GhostLinearSolver;
  class GhostMatrix;
  class GhostVector;
  class StringVector;
}
namespace perulangan
{
  class IterativeSolverVisitorInterface : public alat::InterfaceBase
  {
protected:
    std::string getInterfaceName() const;

public:
    ~IterativeSolverVisitorInterface();
    IterativeSolverVisitorInterface();
    IterativeSolverVisitorInterface( const IterativeSolverVisitorInterface& iterativesolvervisitorinterface);
    IterativeSolverVisitorInterface& operator=( const IterativeSolverVisitorInterface& iterativesolvervisitorinterface);
    std::string getClassName() const;
    perulangan::IterativeSolverVisitorInterface* clone() const;

    // virtual void basicInit(const alat::ParameterFile* parameterfile, std::string blockname);
    virtual std::ostream& printLoopInformation(std::ostream& os) const;
    virtual std::string getVectorType() const;
    // virtual int getVectorLevel() const=0;
    virtual void newVector(alat::GhostVector* u) = 0;

    virtual void vectorEqual(alat::GhostVector& r, const alat::GhostVector& f) const;
    virtual void vectorZero(alat::GhostVector& v) const;
    virtual void vectorAdd(alat::GhostVector& p, double d, const alat::GhostVector& q) const;
    virtual void vectorScale(alat::GhostVector& r, double d) const;
    virtual double vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const;
    virtual double vectorNorm(const alat::GhostVector& r) const;

    virtual void residual(const alat::GhostMatrix& A, alat::GhostVector& r, const alat::GhostVector& u, const alat::GhostVector& f) const;
    virtual void matrixVectorProduct(const alat::GhostMatrix& A, alat::GhostVector& r, const alat::GhostVector& u, double d) const;

    virtual void postProcess(alat::GhostVector& u) const;

    // virtual const alat::armaivec& getDomainsPermutation(int iteration) const;
    // virtual void solveOnDomain(int idomain, const alat::GhostLinearSolver& linearsolverdomain, const alat::GhostMatrix& ghostmatrix, alat::GhostVector& u, const alat::GhostVector& f) const;
    // virtual void vectorEqualOnDomain(int idomain, alat::GhostVector& u, const alat::GhostVector& f) const;
    // virtual void matrixVectorProductCoupling(int i, const alat::GhostMatrix& ghostmatrix, alat::GhostVector& u, const alat::GhostVector& f, double d) const;
    // virtual void smoothInterface(int idomain, alat::GhostVector& u) const;
    // virtual void smoothInterfaceOnLevel(int level, alat::GhostVector& u) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

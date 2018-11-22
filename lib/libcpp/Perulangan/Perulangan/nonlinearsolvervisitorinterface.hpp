#ifndef __Perulangan_NonlinearSolverVisitorInterface_h
#define __Perulangan_NonlinearSolverVisitorInterface_h

#include  "Alat/interfacebase.hpp"
#include  "enums.hpp"

/*--------------------------------------------------------------------------*/

namespace alat
{
  class GhostLinearSolver;
  class GhostMatrix;
  class GhostVector;
  class ParameterFile;
}
namespace perulangan
{
  class NonlinearSolverVisitorInterface : public alat::InterfaceBase
  {
protected:
    std::string getInterfaceName() const;

public:
    ~NonlinearSolverVisitorInterface();
    NonlinearSolverVisitorInterface();
    NonlinearSolverVisitorInterface( const NonlinearSolverVisitorInterface& nonlinearsolvervisitorinterface);
    NonlinearSolverVisitorInterface& operator=( const NonlinearSolverVisitorInterface& nonlinearsolvervisitorinterface);
    std::string getClassName() const;
    NonlinearSolverVisitorInterface* clone() const;

    virtual void init()=0;
    virtual std::ostream& printLoopInformation(std::ostream& os) const;
    virtual std::string getVectorType() const=0;
    virtual void newVector(alat::GhostVector* u) = 0;
    virtual void residual(perulanganEnums::residualstatus& status, alat::GhostVector& gr, const alat::GhostVector& gu, const alat::GhostVector& gf) const=0;
    virtual int solveLinear(perulanganEnums::iterationstatus& status, const alat::GhostLinearSolver& B, const alat::GhostMatrix& A, alat::GhostVector& du, const alat::GhostVector& r) const=0;
    virtual void setLinearTolerance(double rtol, double gtol, alat::GhostLinearSolver& B)=0;

    virtual void vectorZero(alat::GhostVector& gu) const=0;
    virtual void vectorAdd(alat::GhostVector& gu, double s, const alat::GhostVector& gv) const=0;
    virtual double vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const=0;
    virtual void vectorEqual(alat::GhostVector& gu, const alat::GhostVector& gv) const=0;
    virtual std::ostream& vectorWrite(std::ostream& os, const alat::GhostVector& r) const = 0;
    virtual void constructMatrixAndLinearSolvers(perulanganEnums::matrixstatus& status, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u)=0;

    // used for Wolfee line-search
    virtual void computeLinearization(perulanganEnums::residualstatus& status, alat::GhostVector& h, const alat::GhostVector& u, const alat::GhostVector& du) const;
    virtual void setLavrentievParameter(double parameter) const;
    virtual double computeNormSquaredLavrientiev(perulanganEnums::residualstatus& status, const alat::GhostVector& u, const alat::GhostVector& du) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

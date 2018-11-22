#ifndef __Solvers_NonlinearSolverVisitor_hpp
#define __Solvers_NonlinearSolverVisitor_hpp

#include  "Perulangan/enums.hpp"
#include  "Perulangan/nonlinearsolvervisitorinterface.hpp"
#include  "Solvers/solverinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class GhostMatrix;
  class GhostLinearSolver;
  class GhostVector;
}
namespace solvers
{
  class NonlinearSolverVisitor : public perulangan::NonlinearSolverVisitorInterface
  {
  protected:
    SolverInterface* _solver;

  public:
    ~NonlinearSolverVisitor();
    NonlinearSolverVisitor(SolverInterface* solver);
    NonlinearSolverVisitor( const NonlinearSolverVisitor& nonlinearsolvervisitor);
    NonlinearSolverVisitor& operator=( const NonlinearSolverVisitor& nonlinearsolvervisitor);
    std::string getClassName() const;

    void init();
    std::string getVectorType() const;
    void newVector(alat::GhostVector* u);
    void residual(perulanganEnums::residualstatus& status, alat::GhostVector& gr, const alat::GhostVector& gu, const alat::GhostVector& gf) const;
    void setLinearTolerance(double rtol, double gtol, alat::GhostLinearSolver& B);
    void vectorZero(alat::GhostVector& gu) const;
    void vectorAdd(alat::GhostVector& gu, double s, const alat::GhostVector& gv) const;
    double vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const;
    void vectorEqual(alat::GhostVector& gu, const alat::GhostVector& gv) const;
    std::ostream& vectorWrite(std::ostream& os, const alat::GhostVector& r) const;
    int solveLinear(perulanganEnums::iterationstatus& status, const alat::GhostLinearSolver& B, const alat::GhostMatrix& A, alat::GhostVector& du, const alat::GhostVector& r) const;
    void constructMatrixAndLinearSolvers(perulanganEnums::matrixstatus& status, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u);
  };
}

/*--------------------------------------------------------------------------*/
#endif

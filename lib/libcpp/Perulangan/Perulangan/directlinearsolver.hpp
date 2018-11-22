#ifndef __Perulangan_DirectLinearSolver_hpp
#define __Perulangan_DirectLinearSolver_hpp

#include  "Alat/umfmatrix.hpp"
#include  "Alat/matrixallvariables.hpp"
#include  "Perulangan/linearsolverinterface.hpp"
#include  "Solvers/meshunitwithdatainterface.hpp"

/*--------------------------------------------------------------------------*/
namespace perulangan
{
  class DirectLinearSolver : public perulangan::LinearSolverInterface
  {
  protected:
    alat::UmfMatrix _umfmatrix;
    const alat::MatrixAllVariables* _matrixallvariables;
    const solvers::MeshUnitWithDataInterface* _solver;
    alat::MatrixInOne _matrixinone;

  public:
    ~DirectLinearSolver();
    DirectLinearSolver(const alat::MatrixAllVariables* matrix, const solvers::MeshUnitWithDataInterface* solver);
    DirectLinearSolver( const DirectLinearSolver& directlinearsolver);
    DirectLinearSolver& operator=( const DirectLinearSolver& directlinearsolver);
    std::string getClassName() const;

    bool hasIterationInfo() const;
    void reInit();
    void compute();
    void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
    void setTolerance(double rtol, double gtol) {}
  };
}

/*--------------------------------------------------------------------------*/
#endif

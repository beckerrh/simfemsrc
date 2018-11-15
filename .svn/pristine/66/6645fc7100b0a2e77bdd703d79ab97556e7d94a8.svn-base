#ifndef __Solvers_SolverInterface_hpp
#define __Solvers_SolverInterface_hpp

#include  "Alat/interfacebase.hpp"
#include  "Alat/map.hpp"
#include  "Perulangan/enums.hpp"
#include  "Solvers/enums.hpp"
#include  "Solvers/options.hpp"
#include  <memory>
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class GhostMatrix;
  class GhostVector;
  class GhostLinearSolver;
  class MatrixOneVariableInterface;
  class VectorOneVariableInterface;
}
namespace perulangan
{
  class NonlinearSolverInterface;
}
namespace mesh
{
  class MeshInterface;
  class MeshUnitInterface;
}
namespace solvers
{
  class ApplicationInterface;
  class ModelInterface;
  class FemInterface;
  class PdePartInterface;
  class MeshUnitWithDataInterface;
  class Variable;
}

namespace solvers
{
  struct Parameters
  {
    alat::Map<std::string, std::string> strings;
    alat::Map<std::string, int> ints;
    alat::Map<std::string, double> doubles;
    alat::Map<std::string, bool> bools;
  };

  struct varinfo
  {
    const int ncomp;
    const std::string name;
    solverEnums::fem::femtype femname;
    inline explicit varinfo(const std::string name_in, const int ncomp_in, solverEnums::fem::femtype femname_in) : ncomp(ncomp_in), name(name_in), femname(femname_in){}
    bool operator<(const varinfo& v) const;
  };
  struct solverdata
  {
    double timestepratio;
    std::string newton;
    solverdata();
  };

  typedef alat::Map<std::string, arma::vec> ErrorsMap;

  class SolverInterface : public virtual alat::InterfaceBase
  {
  protected:
    virtual std::unique_ptr<solvers::MeshUnitWithDataInterface> newMeshUnitWithDataPlain(const mesh::MeshUnitInterface* mesh) const=0;
    virtual std::unique_ptr<solvers::MeshUnitWithDataInterface> newMeshUnitWithDataBoundary(int color, const mesh::MeshUnitInterface* mesh) const=0;
    virtual std::unique_ptr<solvers::MeshUnitWithDataInterface> newMeshUnitWithDataInterface(int color, const mesh::MeshUnitInterface* mesh) const=0;
    virtual std::unique_ptr<perulangan::NonlinearSolverInterface> newNonlinearSolver() const=0;
    virtual void defineGhosts()=0;

  public:
    ~SolverInterface();
    SolverInterface();
    SolverInterface( const SolverInterface& solverinterface);
    SolverInterface& operator=( const SolverInterface& solverinterface);
    std::string getClassName() const;
    SolverInterface* clone() const;
    virtual const solvers::Parameters& getParameters() const=0;
    virtual std::unique_ptr<alat::MatrixOneVariableInterface> newMatrix() const=0;
    virtual std::unique_ptr<alat::VectorOneVariableInterface> newVector() const=0;
    // virtual std::unique_ptr<solvers::PdePartInterface> newPdePart(const solvers::Variable& var) const=0;
    virtual std::unique_ptr<solvers::FemInterface> newFem(const std::string& varname, solverEnums::fem::femtype fem) const=0;
    virtual std::unique_ptr<solvers::ApplicationInterface> newApplication() const=0;
    virtual std::unique_ptr<solvers::ModelInterface> newModel() const=0;
    virtual void defineVariablesAndPdeParts()=0;
    virtual void enrolVector(alat::GhostVector ghost)=0;
    virtual void enrolMatrix(alat::GhostMatrix ghost)=0;
    virtual void residual(perulanganEnums::residualstatus& status, alat::GhostVector& gr, const alat::GhostVector& gu, const alat::GhostVector& gf) const=0;
    virtual void constructMatrixAndLinearSolvers(perulanganEnums::matrixstatus& status, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u)=0;
    virtual int solveLinear(perulanganEnums::iterationstatus& status, const alat::GhostLinearSolver& B, const alat::GhostMatrix& A, alat::GhostVector& du, const alat::GhostVector& r) const=0;
    virtual solvers::ErrorsMap computeErrors(std::string varname, solver_options::errors::opts opts) const=0;
    virtual void initSolution(alat::GhostVector gu) const=0;
    virtual void initData() const=0;
    virtual void saveSolution(const alat::GhostVector ghost, int it=-1) const=0;
    virtual void saveData(int it=0) const=0;

    virtual perulangan::NonlinearSolverInterface* getNonlinearSolver()=0;
    virtual const perulangan::NonlinearSolverInterface* getNonlinearSolver() const=0;
    virtual void setLinearTolerance(double rtol, double gtol, alat::GhostLinearSolver& B)=0;
    virtual void vectorZero(alat::GhostVector& gu) const=0;
    virtual void vectorAdd(alat::GhostVector& gu, double s, const alat::GhostVector& gv) const=0;
    virtual void vectorEqual(alat::GhostVector& gu, const alat::GhostVector& gv) const=0;
    virtual double vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const=0;
  };
}

/*--------------------------------------------------------------------------*/
#endif

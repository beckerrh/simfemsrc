#ifndef __Solvers_Solver_h
#define __Solvers_Solver_h

#include  "Alat/matrixallvariables.hpp"
#include  "Alat/vectorallvariables.hpp"
#include  "Alat/matrix.hpp"
#include  "Alat/stringvector.hpp"
#include  "Alat/vector.hpp"
#include  "Alat/chronometer.hpp"
#include  "Alat/list.hpp"
#include  "Solvers/meshunitwithdatainterface.hpp"
#include  "Solvers/solverinterface.hpp"
#include  "Solvers/variable.hpp"
#include  "Mesh/meshinterface.hpp"
#include  "Mesh/timemesh.hpp"
#include  "Perulangan/nonlinearsolverinterface.hpp"
#include  <limits>

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class MeasureOfCell;
  class Normals;
}
namespace solvers
{
  class Solver : public virtual solvers::SolverInterface
  {
  public:
    typedef alat::Map<int, std::shared_ptr<solvers::MeshUnitWithDataInterface> >  MeshUnitWithDataMap;

  private:
    int _nvars;
    alat::List<varinfo> _varinfoplain, _varinfodataplain;
    alat::Map<int, varinfo> _varinfobdry, _varinfoitfc;
    alat::Map<std::string, int> _var2ncomp;
    alat::IntSet _bdryunitscolors;

  protected:
    int _printlevel;
    solver_options::opts _opts;
    solver_options::output_manager _output_manager;
    solvers::Parameters _parameters;
    int _partion_id;
    mutable bool _initcalled, _meshsaved;
    mutable alat::Chronometer _chronometer;

    mutable std::shared_ptr<mesh::MeshInterface> _mesh;
    mesh::TimeMesh _timemesh;
    solvers::solverdata _solverdata;

    std::shared_ptr<solvers::MeshUnitWithDataInterface> _plainmeshunitwithdata;
    MeshUnitWithDataMap _boundarymeshunitswithdata;
    MeshUnitWithDataMap _interfacemeshunitswithdata;

    std::shared_ptr<perulangan::NonlinearSolverInterface> _nonlinearsolver;

    std::unique_ptr<solvers::MeshUnitWithDataInterface> newMeshUnitWithDataPlain(const mesh::MeshUnitInterface* mesh) const;
    std::unique_ptr<solvers::MeshUnitWithDataInterface> newMeshUnitWithDataBoundary(int color, const mesh::MeshUnitInterface* mesh) const;
    std::unique_ptr<solvers::MeshUnitWithDataInterface> newMeshUnitWithDataInterface(int color, const mesh::MeshUnitInterface* mesh) const;
    std::unique_ptr<perulangan::NonlinearSolverInterface> newNonlinearSolver() const;

    void initDataDir();
    void addDataVariablePlain(std::string name, int ncomp, solverEnums::fem::femtype fem);
    void addVariablePlain(std::string name, int ncomp, solverEnums::fem::femtype fem);
    void addPdePartPlain(std::string name, std::shared_ptr<solvers::PdePartInterface> pdepart);
    void setTimeMeshData(double time, double dt);
    void defineGhosts();

  public:
    ~Solver();
    Solver( const Solver& solver);
    Solver(std::shared_ptr<mesh::MeshInterface> mesh=nullptr, solver_options::opts opts=solver_options::none, int printlevel=0);
    Solver& operator=( const Solver& solver);
    std::string getClassName() const;
    Solver* clone() const;

    const solvers::Parameters& getParameters() const {return _parameters;}

    void setMesh(std::shared_ptr<mesh::MeshInterface> mesh);
    void setParameter(std::string name, std::string value);
    void setParameter(std::string name, int value);
    void setParameter(std::string name, double value);
    void setParameter(std::string name, bool value);

    std::unique_ptr<alat::MatrixOneVariableInterface> newMatrix() const;
    std::unique_ptr<alat::VectorOneVariableInterface> newVector() const;
    // std::unique_ptr<solvers::PdePartInterface> newPdePart(const solvers::Variable& var) const;
    std::unique_ptr<solvers::FemInterface> newFem(const std::string& varname, solverEnums::fem::femtype fem) const;
    std::unique_ptr<solvers::ApplicationInterface> newApplication() const;
    std::unique_ptr<solvers::ModelInterface> newModel() const;
    void defineVariablesAndPdeParts();

    std::string getInfo() const;
    std::string getMeshInfo() const;
    std::string& setOutputOptions(solver_options::output_manager_data opt);
    void setTimeMesh(const mesh::TimeMeshData timemeshdata);
    void setSolverData(const solvers::solverdata solverdata);

    void loadMesh(std::string meshtype, std::string meshname);
    void loadSolution(std::string filename);
    void saveSolution(const alat::GhostVector ghost, int it=-1) const;
    void saveData(int it=0) const;
    void writeXdmf() const;
    void writeDataXdmf(int it) const;
    void addBoundaryUnits(int first=std::numeric_limits<int>::max(), int last=std::numeric_limits<int>::lowest());
    void init();

    perulangan::NonlinearSolverInterface* getNonlinearSolver();
    const perulangan::NonlinearSolverInterface* getNonlinearSolver() const;
    void enrolVector(alat::GhostVector ghost);
    void enrolMatrix(alat::GhostMatrix ghost);
    void residual(perulanganEnums::residualstatus& status, alat::GhostVector& gr, const alat::GhostVector& gu, const alat::GhostVector& gf) const;
    void constructMatrixAndLinearSolvers(perulanganEnums::matrixstatus& status, alat::GhostLinearSolver& B, alat::GhostMatrix& A, alat::GhostVector& u);
    int solveLinear(perulanganEnums::iterationstatus& status, const alat::GhostLinearSolver& B, const alat::GhostMatrix& A, alat::GhostVector& du, const alat::GhostVector& r) const;
    void setLinearTolerance(double rtol, double gtol, alat::GhostLinearSolver& B);
    void vectorZero(alat::GhostVector& gu) const;
    void vectorAdd(alat::GhostVector& gu, double s, const alat::GhostVector& gv) const;
    void vectorEqual(alat::GhostVector& gu, const alat::GhostVector& gv) const;
    double vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const;
    void initSolution(alat::GhostVector gu) const;
    void initData() const;
    void computeRhs(alat::GhostVector gf, const alat::GhostVector gu) const;
    solvers::ErrorsMap computeErrors(std::string varname, solver_options::errors::opts opts) const;
    void run();
    int macroStep(double time_old, double time, double dtfirst, alat::GhostMatrix gA, alat::GhostLinearSolver gB, alat::GhostVector gu, alat::GhostVector gf, perulangan::NewtonOutputData& newtonoutputdata);
 };
 std::ostream &operator<<(std::ostream &os, const Solver& solver);
}

/*--------------------------------------------------------------------------*/
#endif

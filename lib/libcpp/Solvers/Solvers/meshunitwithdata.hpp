#ifndef __Solvers_MeshUnitWithData_h
#define __Solvers_MeshUnitWithData_h

#include  "Alat/agencyvector.hpp"
#include  "Alat/agencymatrix.hpp"
#include  "Alat/agencylinearsolver.hpp"
#include  "Alat/interfacebase.hpp"
#include  "Alat/matrixallvariables.hpp"
#include  "Alat/vectorallvariables.hpp"
#include  "Alat/matrix.hpp"
#include  "Alat/stringvector.hpp"
#include  "Alat/vector.hpp"
#include  "Alat/chronometer.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  "Solvers/enums.hpp"
#include  "Solvers/pdepartinterface.hpp"
#include  "Solvers/variable.hpp"
#include  "Solvers/meshunitwithdatainterface.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class MeasureOfCell;
  class Normals;
}
namespace solvers
{
  class MeshUnitWithData : public virtual MeshUnitWithDataInterface
  {
  protected:
    mesh::MeshUnitInterface* _mesh;
    const SolverInterface* _solver;
    PdePartsMap _pdeparts;
    bool _initcalled;

    alat::Vector<Variable> _vars, _varsdata;
    alat::StringVector _varnames, _varnamesdata;
    alat::Map<std::string, int> _var2index, _var2indexdata;

    mutable FemMap _fems, _femsex, _femsdata;

    std::shared_ptr<ApplicationInterface> _application;
    std::shared_ptr<ModelInterface> _model;

    mutable alat::AgencyVector _vectoragency;
    mutable alat::AgencyMatrix _matrixagency;
    mutable alat::AgencyLinearSolver _linearsolveragency;

    std::shared_ptr<MeshInfo> _meshinfo;
    mutable PdePartData _pdepartdata;

    mutable alat::armaivec _ncomps, _nlocals;
    void preparePdePartsData(int iK, const alat::armaivec& ivars) const;
    void preparePdePartsData(int iKin, int iKex, const alat::armaivec& ivars) const;
    void prepareCellMatrix(int iK, const alat::armaivec& ivars) const;
    void prepareCellMatrix(int iKin, int iKex, const alat::armaivec& ivars) const;

  public:
    ~MeshUnitWithData();
    MeshUnitWithData();
    MeshUnitWithData( const MeshUnitWithData& meshunitwithdata);
    MeshUnitWithData& operator=( const MeshUnitWithData& meshunitwithdata);
    std::string getClassName() const;
    std::string getInfo() const;
    const alat::Vector<Variable>& getVars() const;
    const alat::Vector<Variable>& getVarsData() const;
    const PdePartsMap& getPdeParts() const;
    PdePartsMap& getPdeParts();
    void setTimeMeshData(double time, double dt);

    const solvers::ApplicationInterface* getApplication() const;
    const solvers::ModelInterface* getModel() const;
    solvers::FemMap* getFemMap() const;
    solvers::FemMap* getFemExMap() const;
    const mesh::MeshUnitInterface* getMesh() const;

    void loadSolution(alat::GhostVector ghost, std::string filename);
    void saveSolution(const alat::GhostVector ghost, std::string filename) const;
    void saveData(std::string filename) const;
    void enrolVector(alat::GhostVector ghost);
    void enrolMatrix(alat::GhostMatrix ghost);
    void enrolLinearSolver(alat::GhostLinearSolver ghost);
    const solvers::FemInterface* getFemData(std::string name) const;
    alat::VectorAllVariables* getVector(alat::GhostVector ghost) const;
    alat::VectorOneVariable* getDataVector(std::string name) const;
    alat::MatrixAllVariables* getMatrix(alat::GhostMatrix ghost) const;
    perulangan::LinearSolverInterface* getLinearSolver(alat::GhostLinearSolver ghost) const;
    void addDataVariable(std::string name, int ncomp, solverEnums::fem::femtype fem);
    void addVariable(std::string name, int ncomp, solverEnums::fem::femtype fem);
    void addPdePart(std::string name, std::shared_ptr<solvers::PdePartInterface> pdepart);
    void initMeshAndApplication(mesh::MeshUnitInterface* mesh, const solvers::SolverInterface* solver);
    void initFemAndMemoryAndDataAndPdeParts();
    void initSolution(alat::GhostVector gu) const;
    void initData() const;
    void initDataVector(alat::VectorOneVariableInterface* u, int ivar, std::string varname) const;

    void residual(perulanganEnums::residualstatus& status, alat::GhostVector& gr, const alat::GhostVector& gu, const alat::GhostVector& gf) const;

    void vectorZero(alat::GhostVector& gu) const;
    void vectorAdd(alat::GhostVector& gu, double s, const alat::GhostVector& gv) const;
    void vectorEqual(alat::GhostVector& gu, const alat::GhostVector& gv) const;
    double vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const;

    void computeRhs(alat::GhostVector gf, const alat::GhostVector gu) const;
    void computeJacobian(perulanganEnums::matrixstatus& status, alat::GhostMatrix gA, const alat::GhostVector gu) const;
    void computeLinearSolver(perulanganEnums::matrixstatus& status, alat::GhostLinearSolver gB, alat::GhostMatrix gA, const alat::GhostVector gu) const;
    int solveLinear(perulanganEnums::iterationstatus& status, const alat::GhostLinearSolver& B, const alat::GhostMatrix& A, alat::GhostVector& du, const alat::GhostVector& r) const;
    void computeErrors(solvers::ErrorsMap& errormaps, std::string varname, const alat::GhostVector gu) const;
 };
}

/*--------------------------------------------------------------------------*/
#endif

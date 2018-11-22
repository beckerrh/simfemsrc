#ifndef __Solvers_MeshUnitWithDataInterface_hpp
#define __Solvers_MeshUnitWithDataInterface_hpp

#include  "Alat/interfacebase.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  "Perulangan/enums.hpp"
#include  "Solvers/solverinterface.hpp"
#include  "Solvers/applicationinterface.hpp"
#include  "Solvers/modelinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class GhostMatrix;
  class GhostVector;
  class GhostLinearSolver;
  class MatrixOneVariableInterface;
  template <class T>
  class Vector;
  class VectorAllVariables;
  class VectorOneVariable;
}
namespace mesh
{
  class MeshUnitInterface;
}
namespace perulangan
{
  class LinearSolverInterface;
}
namespace solvers
{
  class PdePartInterface;
  class SolverInterface;
  class Variable;
}
namespace solvers
{
  struct MeshInfo
  {
    int nnodes, ncells, nsides, nedges, nnodespercell, nedgespercell, nsidespercell, nnodesperside;
    double  dim;
    const arma::mat&  nodes;
    const alat::armaimat& nodes_of_cells;
    const alat::armaimat& sides_of_cells;
    const alat::armaimat& edges_of_cells;
    const alat::armaimat& nodes_of_sides;
    const alat::armaimat& cells_of_sides;
    const alat::armaimat& nodes_of_edges;
    const alat::armaicube& localnodes_of_edges_of_cell;
    const arma::vec&  measure_of_cells;
    const arma::mat&  normals;
    const arma::fmat& sigma;
    const mesh::MeshUnitInterface::BoundaryInformationMap& bdrymesheunitsmap;

    MeshInfo(const arma::mat& in_nodes, const alat::armaimat& in_nodes_of_cells, const alat::armaimat& in_sides_of_cells, const alat::armaimat& in_cells_of_sides, const alat::armaimat& in_edges_of_cells, const alat::armaimat& in_nodes_of_sides, const alat::armaimat& in_nodes_of_edges, const alat::armaicube& in_localnodes_of_edges_of_cell, const arma::vec& in_measure_of_cells, const arma::mat& in_normals, const arma::fmat& in_sigma, double in_dim, int in_nnodes, int in_ncells, int in_nsides, int in_nedges, int in_nnodespercell, int in_nedgespercell, int in_nsidespercell, int in_nnodesperside, const mesh::MeshUnitInterface::BoundaryInformationMap& in_bdrymesheunitsmap);
  };

  typedef alat::Map<std::string, std::shared_ptr<solvers::PdePartInterface> > PdePartsMap;
  typedef alat::Vector<std::shared_ptr<FemInterface> > FemMap;

  class MeshUnitWithDataInterface : public virtual alat::InterfaceBase
  {
  public:
    ~MeshUnitWithDataInterface();
    MeshUnitWithDataInterface();
    MeshUnitWithDataInterface( const MeshUnitWithDataInterface& meshunitwithdatainterface);
    MeshUnitWithDataInterface& operator=( const MeshUnitWithDataInterface& meshunitwithdatainterface);
    std::string getClassName() const;

    virtual void setTimeMeshData(double time, double dt)=0;
    virtual std::string getInfo() const=0;
    virtual const solvers::FemInterface* getFemData(std::string name) const=0;
    virtual const alat::Vector<Variable>& getVars() const=0;
    virtual const alat::Vector<Variable>& getVarsData() const=0;
    virtual const PdePartsMap& getPdeParts() const=0;
    virtual PdePartsMap& getPdeParts()=0;
    virtual void loadSolution(alat::GhostVector ghost, std::string filename)=0;
    virtual void saveSolution(const alat::GhostVector ghost, std::string filename) const=0;
    virtual void saveData(std::string filename) const=0;
    virtual void enrolVector(alat::GhostVector ghost)=0;
    virtual void enrolMatrix(alat::GhostMatrix ghost)=0;
    virtual void enrolLinearSolver(alat::GhostLinearSolver ghost)=0;
    virtual alat::VectorAllVariables* getVector(alat::GhostVector ghost) const=0;
    virtual alat::VectorOneVariable* getDataVector(std::string name) const=0;
    virtual perulangan::LinearSolverInterface* getLinearSolver(alat::GhostLinearSolver ghost) const=0;
    virtual void addDataVariable(std::string name, int ncomp, solverEnums::fem::femtype fem)=0;
    virtual void addVariable(std::string name, int ncomp, solverEnums::fem::femtype fem)=0;
    virtual void addPdePart(std::string name, std::shared_ptr<solvers::PdePartInterface> pdepart)=0;
    virtual void initMeshAndApplication(mesh::MeshUnitInterface* mesh, const solvers::SolverInterface* solver)=0;
    virtual void initFemAndMemoryAndDataAndPdeParts()=0;
    virtual void initSolution(alat::GhostVector gu) const=0;
    virtual void initData() const=0;
    virtual const solvers::ApplicationInterface* getApplication() const=0;
    virtual const solvers::ModelInterface* getModel() const=0;
    virtual solvers::FemMap* getFemMap() const=0;
    virtual solvers::FemMap* getFemExMap() const=0;
    virtual const mesh::MeshUnitInterface* getMesh() const=0;
    virtual void initDataVector(alat::VectorOneVariableInterface* u, int ivar, std::string varname) const=0;

    virtual void residual(perulanganEnums::residualstatus& status, alat::GhostVector& gr, const alat::GhostVector& gu, const alat::GhostVector& gf) const=0;

    virtual void vectorZero(alat::GhostVector& gu) const=0;
    virtual void vectorAdd(alat::GhostVector& gu, double s, const alat::GhostVector& gv) const=0;
    virtual void vectorEqual(alat::GhostVector& gu, const alat::GhostVector& gv) const=0;
    virtual double vectorDot(const alat::GhostVector& gu, const alat::GhostVector& gv) const=0;
    virtual int solveLinear(perulanganEnums::iterationstatus& status, const alat::GhostLinearSolver& B, const alat::GhostMatrix& A, alat::GhostVector& du, const alat::GhostVector& r) const=0;

    virtual void computeRhs(alat::GhostVector gf, const alat::GhostVector gu) const=0;
    virtual void computeJacobian(perulanganEnums::matrixstatus& status, alat::GhostMatrix gA, const alat::GhostVector gu) const=0;
    virtual void computeLinearSolver(perulanganEnums::matrixstatus& status, alat::GhostLinearSolver gB, alat::GhostMatrix gA, const alat::GhostVector gu) const=0;
    virtual void computeErrors(solvers::ErrorsMap& errormaps, std::string varname, const alat::GhostVector gu) const=0;
  };
}

/*--------------------------------------------------------------------------*/
#endif

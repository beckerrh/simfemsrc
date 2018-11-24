#ifndef ___MeshUnitWithData_hpp
#define ___MeshUnitWithData_hpp

#include  "Solvers/meshunitwithdata.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Alat/cubicspline.hpp"

/*--------------------------------------------------------------------------*/
class CutInterfaceConstructor : public mesh::CutInterfaceConstructor
{
private:
  const alat::VectorOneVariable* _vector;
public:
  CutInterfaceConstructor(const alat::VectorOneVariable* vector);
  std::string getClassName() const;
  const alat::VectorOneVariable* getPhi() const;
};

/*--------------------------------------------------------------------------*/
  class MeshUnitWithData : public solvers::MeshUnitWithData
  {
  private:
    alat::VectorOneVariable _phivector;
    std::string _applicationname;
    alat::CubicSpline _spline;

  protected:
    double _residualcorrectPhi(const solvers::FunctionInterface& phi, const alat::VectorOneVariable& phivector, int& ncutedges) const;
    void _writeVtkCut(std::ofstream& file, const alat::armaivec& nodeintonode, const alat::IntMap& cellintocell, const alat::IntMap& nodetonodein, const alat::armavec& u);

public:
    ~MeshUnitWithData();
    MeshUnitWithData(std::string applicationname);
    MeshUnitWithData( const MeshUnitWithData& meshunitwithdata);
    MeshUnitWithData& operator=( const MeshUnitWithData& meshunitwithdata);
    std::string getClassName() const;

    void initMeshAndApplication(mesh::MeshUnitInterface* mesh, const solvers::SolverInterface* solver);
    void initDataVector(alat::VectorOneVariableInterface* u, int ivar, std::string varname) const;

    void constructPhiByHand();
    void constructPhiBySplines(std::string applicationname);
    void writeVtkCut(std::string varname, std::string filename);
    void writeVtkIIM(std::string varname, std::string filename);
  };

/*--------------------------------------------------------------------------*/
#endif

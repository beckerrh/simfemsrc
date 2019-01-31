#ifndef __Solvers_PdePartInterface_hpp
#define __Solvers_PdePartInterface_hpp

#include  "Alat/interfacebase.hpp"
#include  "Alat/matrix.hpp"
#include  "Alat/strings.hpp"
#include  "Alat/vector.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Solvers/solverinterface.hpp"
#include  "Solvers/feminterface.hpp"
#include  "Solvers/options.hpp"
#include  "Solvers/meshunitwithdatainterface.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class MatrixAllVariables;
  class MatrixOneVariableInterface;
  template <class T>
  class Vector;
  class VectorAllVariables;
}
namespace mesh
{
  class MeshUnitInterface;
}
namespace solvers
{
  struct MeshInfo;
  class Variable;
}
namespace solvers
{
  struct PdePartData
  {
    // typedef arma::field<alat::armaimat> ivec;
    typedef arma::field<alat::armaivec> ivec;
    // typedef arma::field<arma::mat> vec;
    typedef arma::field<alat::armavec> vec;
    // typedef arma::field<alat::armavec> mat;
    // typedef arma::field<alat::armaivec> imat;
    typedef arma::field<arma::mat> mat;
    vec uloc, floc, ulocex, flocex;
    ivec vec_i, vec_iex;
    mat aloc, aloc_inex, aloc_exin, aloc_exex;
    // imat aloc_i, aloc_j, aloc_i_ex, aloc_j_ex;
    void set_size(const alat::armaivec& ncomps, const alat::armaivec& nlocals);
    void set_sizes(const alat::armaivec& ncomps, const alat::armaivec& nlocals);
  };

  typedef alat::Vector< const FemData*> FemDatas;

  class PdePartInterface : public alat::InterfaceBase
  {
  protected:
    const mesh::MeshUnitInterface* _mesh;
    const MeshUnitWithDataInterface* _meshunit;
    const MeshInfo* _meshinfo;
    double _time, _dt;
    const ApplicationInterface* _application;
    const ModelInterface* _model;
    FemMap* _fems, *_femsex;
    mutable FemDatas _femdatas, _femdatasex;
    const FemInterface* _femforintegration;
    alat::StringList _vars;
    alat::armaivec _ivars;
    solver_options::pdepart::opts _opts;
    arma::uvec _cellisbdry;

  public:
    ~PdePartInterface();
    PdePartInterface(alat::StringList vars, solver_options::pdepart::opts opts=solver_options::pdepart::none);
    PdePartInterface( const PdePartInterface& pdepartinterface);
    PdePartInterface& operator=( const PdePartInterface& pdepartinterface);
    std::string getClassName() const;
    std::string getInfo() const;
    const alat::armaivec& getIvars() const;

    virtual const solver_options::pdepart::opts setOptions()=0;
    virtual void initPdePart(const solvers::MeshInfo* meshinfo, const alat::Map<std::string, int>& var2index, const solvers::MeshUnitWithDataInterface* meshunitwithdata, const solvers::Parameters& parameters);
    virtual void setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters);

    bool loopCells() const;
    bool loopBoundary() const;
    bool loopInteriorSides() const;
    virtual bool interiorsidecoupling(int iKin, int iKex) const;
    virtual void rhsCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const;
    virtual void residualCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const;
    virtual void matrixCell(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const;
    virtual void rhsBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const;
    virtual void residualBdry(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const;
    virtual void matrixBdry(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const;
    virtual void rhsBdryCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const;
    virtual void residualBdryCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const;
    virtual void matrixBdryCell(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const;
    virtual void prepareRhsCellBdry(int iK) const;

    virtual void computeResidualInteriorSide(int iS, int iKin, int iKex, solvers::PdePartData::vec& flocin, solvers::PdePartData::vec& flocex, const solvers::PdePartData::vec& ulocin, const solvers::PdePartData::vec& ulocex);
    virtual void computeMatrixInteriorSide(int iS, int iKin, int iKex, solvers::PdePartData::mat& matinin, solvers::PdePartData::mat& matinex, solvers::PdePartData::mat& matexin, solvers::PdePartData::mat& matexex, const solvers::PdePartData::vec& ulocin, const solvers::PdePartData::vec& ulocex)const;
    virtual void computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    virtual void computeRhsBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    virtual void computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    virtual void computeResidualBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    virtual void computeMatrixCell(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const;
    virtual void computeMatrixBdry(int color, int iK, int iS, int iil, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const;
    virtual void setTimeMeshData(double time, double dt);
    virtual void additionCouplings(alat::Matrix<alat::SparsityPatternSoft>& sparsitypatternsoft)const;
    virtual void computeMatrixGlobal(alat::MatrixAllVariables& A, const alat::VectorAllVariables& u)const;
    virtual void computeResidualGlobal(alat::VectorAllVariables& r, const alat::VectorAllVariables& u)const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

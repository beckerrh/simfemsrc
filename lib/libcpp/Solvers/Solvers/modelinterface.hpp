#ifndef __Solvers_ModelInterface_hpp
#define __Solvers_ModelInterface_hpp

#include  "Alat/interfacebase.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  struct FemData;
}
namespace solvers
{
  struct Parameters;
}
namespace solvers
{
  class ModelInterface : public alat::InterfaceBase
  {
  protected:
    const mesh::MeshUnitInterface* _mesh;

  public:
    ~ModelInterface();
    ModelInterface();
    ModelInterface( const ModelInterface& modelinterface);
    ModelInterface& operator=( const ModelInterface& modelinterface);
    std::string getClassName() const;
    virtual std::string getInfo() const=0;

    virtual void initModel(const mesh::MeshUnitInterface* mesh, const solvers::Parameters& parameters);
    // virtual void residualCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems) const=0;
    // virtual void matrixCell(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const=0;
  };
}

/*--------------------------------------------------------------------------*/
#endif

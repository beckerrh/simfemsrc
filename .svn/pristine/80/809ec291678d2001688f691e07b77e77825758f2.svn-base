#ifndef __Solvers_Model_hpp
#define __Solvers_Model_hpp

#include  "modelinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  class Model : public ModelInterface
  {
  protected:
  public:
    ~Model();
    Model();
    Model( const Model& model);
    Model& operator=( const Model& model);
    std::string getClassName() const;
    std::string getInfo() const;

    // void residualCell(solvers::PdePartData::vec& floc, const solvers::FemDatas& fems)const;
    // void matrixCell(solvers::PdePartData::mat& mat, const solvers::FemDatas& fems)const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

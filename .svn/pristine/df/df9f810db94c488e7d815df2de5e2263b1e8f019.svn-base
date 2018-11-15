#ifndef ___Model_hpp
#define ___Model_hpp

#include  "Solvers/modelinterface.hpp"

/*--------------------------------------------------------------------------*/
  class Model : public solvers::ModelInterface
  {
  public:
    ~Model();
    Model();
    Model( const Model& model);
    Model& operator=( const Model& model);
    std::string getClassName() const;
    std::string getInfo() const;

    double a, b, k0, k1;

    void initModel(const mesh::MeshUnitInterface* mesh, const solvers::Parameters& parameters);
  };

/*--------------------------------------------------------------------------*/
#endif

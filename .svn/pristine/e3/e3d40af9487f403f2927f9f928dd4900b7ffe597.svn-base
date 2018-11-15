#ifndef __Solvers_ApplicationInterface_hpp
#define __Solvers_ApplicationInterface_hpp

#include  "Alat/interfacebase.hpp"
#include  "Alat/vector.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  "Solvers/functioninterface.hpp"
#include  "Solvers/modelinterface.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class StringVector;
}
namespace solvers
{

  class ApplicationInterface : public alat::InterfaceBase
  {
  protected:
    const mesh::MeshUnitInterface* _mesh;
    alat::Vector<std::shared_ptr<FunctionInterface> > _exactsolutions;
    const solvers::ModelInterface* _model;
    alat::IntSet _dircolors;

    alat::Vector<std::shared_ptr<InitialConditionInterface> > _initialconditions;
    alat::Vector<std::shared_ptr<RightHandSideInterface> > _righthandsides;
    alat::Vector<std::shared_ptr<DirichletInterface> > _dirichlets;
    alat::Vector<std::shared_ptr<NeumannInterface> > _neumanns;
    alat::Vector<std::shared_ptr<FunctionInterface> > _datafunctions;
    alat::Map<std::string, int> _indexofvar, _indexofdatavar;

    virtual std::unique_ptr<solvers::InitialConditionInterface> newInitialCondition(std::string varname) const;
    virtual std::unique_ptr<solvers::RightHandSideInterface> newRightHandSide(std::string varname) const;
    virtual std::unique_ptr<solvers::DirichletInterface> newDirichlet(std::string varname) const;
    virtual std::unique_ptr<solvers::NeumannInterface> newNeumann(std::string varname) const;
    virtual std::shared_ptr<solvers::FunctionInterface> newDataFunction(std::string varname) const;
    virtual std::shared_ptr<solvers::FunctionInterface> newExactSolution(std::string varname) const;

  public:
    ~ApplicationInterface();
    ApplicationInterface();
    ApplicationInterface( const ApplicationInterface& applicationinterface);
    ApplicationInterface& operator=( const ApplicationInterface& applicationinterface);
    std::string getClassName() const;
    virtual std::string getInfo() const;
    void setModel(const solvers::ModelInterface* model);

    virtual void setExactSolutions(solvers::FunctionInterface& exactsolution) {}
    virtual void initApplication(const mesh::MeshUnitInterface* mesh, const alat::StringVector& varnames, const alat::StringVector& varnamesdata, const alat::armaivec& ncomps, const alat::armaivec& ncompsdata);

    const InitialConditionInterface& getInitialCondition(int ivar)const;
    const RightHandSideInterface& getRightHandSide(int ivar)const;
    const DirichletInterface& getDirichlet(int ivar)const;
    const NeumannInterface& getNeumann(int ivar)const;
    const FunctionInterface& getDataFunction(int ivar)const;
    const solvers::FunctionInterface& getExactSolution(int ivar)const;
    solvers::FunctionInterface& getExactSolution(int ivar);

    const InitialConditionInterface& getInitialCondition(std::string varname)const;
    const RightHandSideInterface& getRightHandSide(std::string varname)const;
    const DirichletInterface& getDirichlet(std::string varname)const;
    const NeumannInterface& getNeumann(std::string varname)const;
    const FunctionInterface& getDataFunction(std::string varname)const;
    const solvers::FunctionInterface& getExactSolution(std::string varname)const;
    solvers::FunctionInterface& getExactSolution(std::string varname);

    bool hasInitialCondition(int ivar)const;
    bool hasRightHandSide(int ivar)const;
    bool hasDirichlet(int ivar)const;
    bool hasDataFunction(int ivar)const;
    bool hasExactSolution(int ivar)const;
    bool hasExactSolution(std::string varname)const;

    virtual bool isStrongDirichlet(int color)const=0;
    virtual const alat::IntSet& getStrongDirichletColor() const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

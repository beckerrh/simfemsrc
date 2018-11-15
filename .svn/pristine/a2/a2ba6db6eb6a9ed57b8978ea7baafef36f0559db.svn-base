#ifndef __Solvers_Variable_hpp
#define __Solvers_Variable_hpp

#include  "Alat/interfacebase.hpp"
#include  "Solvers/enums.hpp"

/*--------------------------------------------------------------------------*/
namespace mesh
{
  class MeshUnitInterface;
}
namespace solvers
{
  class Variable : public virtual alat::InterfaceBase
  {
  protected:
    std::string _name;
    solverEnums::fem::femtype _fem;
    int _ncomp;

  public:
    ~Variable();
    Variable();
    Variable(std::string name, int ncomp, solverEnums::fem::femtype fem);
    Variable( const Variable& variable);
    Variable& operator=( const Variable& variable);
    std::string getClassName() const;
    Variable* clone() const;

    std::string getName() const;
    int getNcomp() const;
    solverEnums::fem::femtype getFem() const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

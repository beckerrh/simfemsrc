#ifndef __Solvers_Timestepping_hpp
#define __Solvers_Timestepping_hpp

#include  "Alat/interfacebase.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  class Timestepping : public alat::InterfaceBase
  {
  private:
  protected:
  public:
    ~Timestepping();
    Timestepping();
    Timestepping( const Timestepping& timestepping);
    Timestepping& operator=( const Timestepping& timestepping);
    std::string getClassName() const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

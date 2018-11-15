#include  "Solvers/timestepping.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
Timestepping::~Timestepping() {}
Timestepping::Timestepping(): alat::InterfaceBase(){}
Timestepping::Timestepping( const Timestepping& timestepping): alat::InterfaceBase(timestepping)
{
  (*this).operator=(timestepping);
}
Timestepping& Timestepping::operator=( const Timestepping& timestepping)
{
  assert(0);
  alat::InterfaceBase::operator=(timestepping);
  return *this;
}
std::string Timestepping::getClassName() const 
{
  return "Timestepping";
}
/*--------------------------------------------------------------------------*/

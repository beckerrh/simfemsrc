#include  "Solvers/iodata.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
IoData::~IoData() {}
IoData::IoData(): alat::InterfaceBase(){}
IoData::IoData( const IoData& iodata): alat::InterfaceBase(iodata)
{
  (*this).operator=(iodata);
}
IoData& IoData::operator=( const IoData& iodata)
{
  assert(0);
  alat::InterfaceBase::operator=(iodata);
  return *this;
}
std::string IoData::getClassName() const 
{
  return "IoData";
}
/*--------------------------------------------------------------------------*/

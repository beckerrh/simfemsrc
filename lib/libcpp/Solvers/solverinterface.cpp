#include  "Solvers/solverinterface.hpp"
#include  <cassert>

using namespace solvers;

/*--------------------------------------------------------------------------*/
bool varinfo::operator<(const varinfo& v) const
{
  if(ncomp<v.ncomp) return true;
  if(ncomp>v.ncomp) return false;
  if(name<v.name) return true;
  if(name>v.name) return false;
  if(femname<v.femname) return true;
  if(femname>v.femname) return false;
  assert(0); return true;
}
/*--------------------------------------------------------------------------*/
  solverdata::solverdata()
  {
    // nlmaxiter=10;
    // nlrtol=1e-6;
    timestepratio=1.0;
    // rhomatrix=0.025;
    newton="newton";
  }

/*--------------------------------------------------------------------------*/
SolverInterface::~SolverInterface() {}
SolverInterface::SolverInterface(): alat::InterfaceBase(){}
SolverInterface::SolverInterface( const SolverInterface& solverinterface): alat::InterfaceBase(solverinterface)
{
  (*this).operator=(solverinterface);
}
SolverInterface& SolverInterface::operator=( const SolverInterface& solverinterface)
{
  assert(0);
  alat::InterfaceBase::operator=(solverinterface);
  return *this;
}
std::string SolverInterface::getClassName() const
{
  return "SolverInterface";
}
SolverInterface* SolverInterface::clone() const
{
	assert(0);
	return NULL;
  // return new SolverInterface(*this);
}
/*--------------------------------------------------------------------------*/

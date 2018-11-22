#include  "Perulangan/iterationinfo.hpp"
#include  "Perulangan/linearsolverinterface.hpp"
#include  <cassert>
#include  <iostream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
LinearSolverInterface::~LinearSolverInterface(){}
LinearSolverInterface::LinearSolverInterface() : alat::InterfaceBase(){}
LinearSolverInterface::LinearSolverInterface( const LinearSolverInterface& linearsolverinterface) : alat::InterfaceBase(linearsolverinterface){}
LinearSolverInterface& LinearSolverInterface::operator=( const LinearSolverInterface& linearsolverinterface)
{
  InterfaceBase::operator=(linearsolverinterface);
  assert(0);
  return *this;
}
std::ostream& LinearSolverInterface::printLoopInformation(std::ostream& os) const
{
  os << "\t" << getClassName() << "\n";
  return os;
}
std::string LinearSolverInterface::getInterfaceName() const
{
  return "LinearSolverInterface";
}

/*--------------------------------------------------------------------------*/
int LinearSolverInterface::getNVectors() const
{
  return 0;
}

/*--------------------------------------------------------------------------*/
std::ostream& LinearSolverInterface::write(std::ostream& os) const
{
  _notWritten("write");
  return os;
}

/*--------------------------------------------------------------------------*/
// void LinearSolverInterface::basicInit(const alat::ParameterFile* parameterfile, std::string blockname)
// {
//   _notWritten("basicInit");
// }

// void LinearSolverInterface::reInit(){_notWritten("reInit");}
// void LinearSolverInterface::compute(){_notWritten("compute");}
void LinearSolverInterface::addUpdate(perulanganEnums::iterationstatus& status, const alat::GhostVector& w, alat::GhostVector& u) const{_notWritten("addUpdate");}
const IterationInfo* LinearSolverInterface::getIterationInfo() const {_notWritten("getIterationInfo"); return NULL;}
IterationInfo* LinearSolverInterface::getIterationInfo() {_notWritten("getIterationInfo"); return NULL;}

/*--------------------------------------------------------------------------*/
void LinearSolverInterface::restart()
{
  // _notWritten("restart");
}

/*--------------------------------------------------------------------------*/
void LinearSolverInterface::setTolerance(double rtol, double gtol)
{
  IterationInfoData& data = getIterationInfo()->getIterationInfoData();
  data.rtol = rtol;
  data.gtol = gtol;
}

int LinearSolverInterface::getNumberOfIterations() const
{
  return getIterationInfo()->getNumberOfIterations();
}

#include  "Perulangan/nonlinearsolverinterface.hpp"
#include  <cassert>
#include  <iostream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
NewtonInputData::NewtonInputData()
{
  printlinesearch = false;
  maxnlinesearch = 10;
  maxiter = 10;
  printstep = 1;
  gtol = 1e-12;
  rtol = 1e-5;
  omegalinesearch = 0.7;
  monotonyfactor = 1.0;
  lineartoleranceincrease = 0.2;
  rhomatrix = 0.025;
}
NewtonInputData::NewtonInputData(const NewtonInputData& newtoninputdata)
{
  printlinesearch = newtoninputdata.printlinesearch;
  maxnlinesearch = newtoninputdata.maxnlinesearch;
  maxiter = newtoninputdata.maxiter;
  printstep = newtoninputdata.printstep;
  gtol = newtoninputdata.gtol;
  rtol = newtoninputdata.rtol;
  omegalinesearch = newtoninputdata.omegalinesearch;
  monotonyfactor = newtoninputdata.monotonyfactor;
  lineartoleranceincrease = newtoninputdata.lineartoleranceincrease;
  rhomatrix = newtoninputdata.rhomatrix;
}
std::ostream& perulangan::operator<<(std::ostream& os, const NewtonInputData& nd)
{
  os << "lineartoleranceincrease/rhomatrix/maxnlinesearch/omegalinesearch/monotonyfactor " <<nd.lineartoleranceincrease<<"/"<<nd.rhomatrix<<"/"<<nd.maxnlinesearch<<"/"<<nd.omegalinesearch<<"/"<<nd.monotonyfactor<<"\n";
  return os;
}
/*-------------------------------------------------------------*/
NewtonOutputData::NewtonOutputData()
{
  niter = -1;
  niter_linear = -1;
  nredo_matrix = -1;
  residual_status = perulanganEnums::ResidualStatusNone;
  linear_solver_status = perulanganEnums::IterationStatusNone;
  newton_status = perulanganEnums::NewtonStatusNone;
  matrix_status = perulanganEnums::MatrixStatusNone;
}
std::ostream& perulangan::operator<<(std::ostream& os, const NewtonOutputData& n)
{
  os << "(niter/niter_linear/nredo_matrix:) " << n.niter <<"/"<< n.niter_linear <<"/"<< n.nredo_matrix<<"/"<< "\n";
  os << "(residual_status/linear_solver_status/newton_status/matrix_status:) " << perulanganEnums::residualStatusToString(n.residual_status)<<"/"
  << perulanganEnums::iterationStatusToString(n.linear_solver_status) <<"/"<< perulanganEnums::newtonStatusToString(n.newton_status)<<"/" << perulanganEnums::matrixStatusToString(n.matrix_status);
  return os;
}

/*--------------------------------------------------------------------------*/
NonlinearSolverInterface::~NonlinearSolverInterface(){}
NonlinearSolverInterface::NonlinearSolverInterface() : alat::InterfaceBase(){}
NonlinearSolverInterface::NonlinearSolverInterface( const NonlinearSolverInterface& linearsolverinterface) : alat::InterfaceBase(linearsolverinterface){}
NonlinearSolverInterface& NonlinearSolverInterface::operator=( const NonlinearSolverInterface& linearsolverinterface)
{
  InterfaceBase::operator=(linearsolverinterface);
  assert(0);
  return *this;
}
std::string NonlinearSolverInterface::getInterfaceName() const
{
  return "NonlinearSolverInterface";
}

/*--------------------------------------------------------------------------*/
std::ostream& NonlinearSolverInterface::printLoopInformation(std::ostream& os) const
{
  os << "\t" << getClassName() << "\n";
  return os;
}
// const IterationInfo& NonlinearSolverInterface::getIterationInfo() const{_notWritten("getIterationInfo");}
// IterationInfo& NonlinearSolverInterface::getIterationInfo(){_notWritten("getIterationInfo");}

/*--------------------------------------------------------------------------*/
int NonlinearSolverInterface::getNVectors() const{_notWritten("getNVectors"); return 0;}
std::ostream& NonlinearSolverInterface::write(std::ostream& os) const{_notWritten("write");return os;}
void NonlinearSolverInterface::reInit(){_notWritten("reInit");}
void NonlinearSolverInterface::addUpdate(perulanganEnums::iterationstatus& status, const alat::GhostVector& w, alat::GhostVector& u) const{_notWritten("addUpdate");}

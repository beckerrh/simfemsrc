#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/iterativesolverwithpreconditioner.hpp"
#include  "Perulangan/preconditionerinterface.hpp"
#include  <cassert>
#include  <sstream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
IterativeSolverWithPreconditioner::~IterativeSolverWithPreconditioner()
{
  if(_deletepreconditioner && _preconditioner)
  {
    delete _preconditioner;
    _preconditioner = NULL;
  }
}

IterativeSolverWithPreconditioner::IterativeSolverWithPreconditioner() : IterativeSolverWithVisitor(), _preconditioner(NULL), _deletepreconditioner(true){}
IterativeSolverWithPreconditioner::IterativeSolverWithPreconditioner( const IterativeSolverWithPreconditioner& iterativesolverwithpreconditioner) : IterativeSolverWithVisitor(iterativesolverwithpreconditioner), _preconditioner(iterativesolverwithpreconditioner._preconditioner), _deletepreconditioner(iterativesolverwithpreconditioner._deletepreconditioner)
{
  assert(0);
}

IterativeSolverWithPreconditioner& IterativeSolverWithPreconditioner::operator=( const IterativeSolverWithPreconditioner& iterativesolverwithpreconditioner)
{
  IterativeSolverWithVisitor::operator=(iterativesolverwithpreconditioner);
  _preconditioner = iterativesolverwithpreconditioner._preconditioner;
  _deletepreconditioner = iterativesolverwithpreconditioner._deletepreconditioner;
  return *this;
}

std::string IterativeSolverWithPreconditioner::getClassName() const
{
  return "IterativeSolverWithPreconditioner";
}

perulangan::IterativeSolverWithPreconditioner*  IterativeSolverWithPreconditioner::clone() const
{
  assert(0);
  // return new IterativeSolverWithPreconditioner(*this);
}

/*--------------------------------------------------------------------------*/
// void IterativeSolverWithPreconditioner::basicInit(const alat::ParameterFile* parameterfile, std::string blockname)
// {
//   // std::cerr << "IterativeSolverWithPreconditioner::basicInit() " << getClassName() << " --> visitor =  " << getVisitor()->getClassName() << " --> preconditioner =  " << getPreconditioner()->getClassName() <<"\n";
//   IterativeSolverWithVisitor::basicInit(parameterfile, blockname);
//   assert(_preconditioner);
//   getPreconditioner()->basicInit( parameterfile, blockname, getVisitor() );
//   std::stringstream ss;
//   ss<<"\t"<< getClassName() << " ("<< getVisitor()->getClassName() << ") : " << getPreconditioner()->getClassName() << " ";
//   getIterationInfo()->setId( ss.str() );
// }

/*--------------------------------------------------------------------------*/
void IterativeSolverWithPreconditioner::reInit()
{
  assert(_basicinitcalled==true);
  getPreconditioner()->reInit();
}

/*--------------------------------------------------------------------------*/
void IterativeSolverWithPreconditioner::compute()
{
  getPreconditioner()->computePreconditioner();
}

/*--------------------------------------------------------------------------*/
std::ostream& IterativeSolverWithPreconditioner::printLoopInformation(std::ostream& os) const
{
  IterativeSolverWithVisitor::printLoopInformation(os);
  os << " ";
  getPreconditioner()->printLoopInformation(os);
  return os;
}

/*--------------------------------------------------------------------------*/
perulangan::PreconditionerInterface*& IterativeSolverWithPreconditioner::newPreconditionerPointer()
{
  assert(_preconditioner == NULL);
  _deletepreconditioner = true;
  return _preconditioner;
}

/*--------------------------------------------------------------------------*/
void IterativeSolverWithPreconditioner::setPreconditionerPointer(perulangan::PreconditionerInterface* preconditioner)
{
  _preconditioner = preconditioner;
  _deletepreconditioner = false;
}

/*--------------------------------------------------------------------------*/
const perulangan::PreconditionerInterface* IterativeSolverWithPreconditioner::getPreconditioner() const
{
  assert(_preconditioner);
  return _preconditioner;
}

/*--------------------------------------------------------------------------*/
perulangan::PreconditionerInterface* IterativeSolverWithPreconditioner::getPreconditioner()
{
  assert(_preconditioner);
  return _preconditioner;
}

/*--------------------------------------------------------------------------*/
void IterativeSolverWithPreconditioner::addUpdate(perulanganEnums::iterationstatus& status, const alat::GhostVector& w, alat::GhostVector& u) const
{
  getVisitor()->vectorAdd( u, 1.0, w );
}

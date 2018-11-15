#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/preconditioner.hpp"
#include  <cassert>
#include  <sstream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
int Preconditioner::_totalids = 0;

Preconditioner::~Preconditioner(){}
Preconditioner::Preconditioner() : PreconditionerInterface(), _visitor(NULL)
{
  _id = _totalids;
  _totalids++;
}
Preconditioner::Preconditioner( const Preconditioner& preconditioner) : PreconditionerInterface(preconditioner), _visitor(preconditioner._visitor)
{
  assert(0);
}
Preconditioner& Preconditioner::operator=( const Preconditioner& preconditioner)
{
  PreconditionerInterface::operator=(preconditioner);
  assert(0);
  return *this;
}
std::string Preconditioner::getClassName() const
{
  return "constructLinearSolver";
}
Preconditioner* Preconditioner::clone() const
{
  assert(0);
//return new Preconditioner(*this);
}

/*--------------------------------------------------------------------------*/
// void Preconditioner::basicInit(const alat::ParameterFile* parameterfile, std::string blockname, perulangan::IterativeSolverVisitorInterface* visitor)
// {
//   // std::cerr << "Preconditioner::basicInit() " << getClassName() << "\n";
//   _visitor = visitor;
//   assert(_visitor);
//   memory();
// }

/*--------------------------------------------------------------------------*/
void Preconditioner::setsmoothtype(std::string smoothtype)
{
  _smoothtype = smoothtype;
}

/*--------------------------------------------------------------------------*/
int Preconditioner::getNVectors() const
{
  return 0;
}

/*--------------------------------------------------------------------------*/
const perulangan::IterativeSolverVisitorInterface* Preconditioner::getVisitor() const
{
  assert(_visitor);
  return _visitor;
}

perulangan::IterativeSolverVisitorInterface* Preconditioner::getVisitor()
{
  assert(_visitor);
  return _visitor;
}

/*--------------------------------------------------------------------------*/
alat::GhostVector& Preconditioner::getMemory(int i) const
{
  if( i >= _memory.size() )
  {
    _error_string("getMemory", "", "too small memory");
  }
  return _memory[i];
}

/*--------------------------------------------------------------------------*/
void Preconditioner::memory()
{
  int nvectors = getNVectors();
  if(nvectors == 0)
  {
    return;
  }
  if(_visitor == NULL)
  {
    _error_string("memory", "no visitor !");
  }
  _memory.set_size( getNVectors() );
  std::string type = getVisitor()->getVectorType();
  // int level = getVisitor()->getVectorLevel();
  // std::cerr << "Preconditioner::memory() in " << getClassName() << " visitor = "<<getVisitor()->getClassName() << " " << getNVectors()  << " of type " << type << "\n";
  // assert(0);
  for(int i = 0; i < _memory.size(); i++)
  {
    std::stringstream ss;
    ss<<getClassName()<<"_"<<_id<<"_"<<getVisitor()->getClassName()<<"_memory_"<<i;
    _memory[i] = alat::GhostVector( ss.str(), type);
    // _memory[i].setVariables( variablesOfVector(i) );
    getVisitor()->newVector(&_memory[i]);
  }
}

#include  "Perulangan/iterativesolvervisitorinterface.hpp"
#include  "Perulangan/iterativesolverwithvisitor.hpp"
#include  <cassert>
#include  <sstream>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
int IterativeSolverWithVisitor::_totalids = 0;

IterativeSolverWithVisitor::~IterativeSolverWithVisitor()
{
  if(_deletevisitor && _visitor)
  {
    delete _visitor;
    _visitor = NULL;
  }
}

IterativeSolverWithVisitor::IterativeSolverWithVisitor() : LinearSolverInterface(), _visitor(NULL), _deletevisitor(true), _basicinitcalled(false), _info()
{
  _id = _totalids;
  _totalids++;
}

IterativeSolverWithVisitor::IterativeSolverWithVisitor( const IterativeSolverWithVisitor& iterativesolverwithvisitor) : LinearSolverInterface(iterativesolverwithvisitor), _deletevisitor(iterativesolverwithvisitor._deletevisitor), _info(iterativesolverwithvisitor._info)
{
  _id = _totalids;
  _totalids++;
  if(_deletevisitor)
  {
    if(not iterativesolverwithvisitor._visitor)
    {
      _visitor = NULL;
    }
    else
    {
      _visitor = iterativesolverwithvisitor._visitor->clone();
    }
  }
  else
  {
    _visitor = iterativesolverwithvisitor._visitor;
  }
}

IterativeSolverWithVisitor& IterativeSolverWithVisitor::operator=( const IterativeSolverWithVisitor& iterativesolverwithvisitor)
{
  LinearSolverInterface::operator=(iterativesolverwithvisitor);
  _visitor = iterativesolverwithvisitor._visitor;
  _deletevisitor = iterativesolverwithvisitor._deletevisitor;
  return *this;
}

std::string IterativeSolverWithVisitor::getClassName() const
{
  return "IterativeSolverWithVisitor";
}

IterativeSolverWithVisitor* IterativeSolverWithVisitor::clone() const
{
  assert(0);
  // return new IterativeSolverWithVisitor(*this);
}

/*--------------------------------------------------------------------------*/
std::ostream& IterativeSolverWithVisitor::printLoopInformation(std::ostream& os) const
{
  os << "\"" << getClassName() << "\" ";
  getVisitor()->printLoopInformation(os);
  os << " info: ";
  getIterationInfo()->printLoopInformation(os);
  os <<" ";
  return os;
}

bool IterativeSolverWithVisitor::hasIterationInfo() const
{
  return true;
}

/*--------------------------------------------------------------------------*/
perulangan::IterativeSolverVisitorInterface*& IterativeSolverWithVisitor::newVisitorPointer()
{
  assert(_visitor == NULL);
  _deletevisitor = true;
  return _visitor;
}

void IterativeSolverWithVisitor::setVisitorPointer(perulangan::IterativeSolverVisitorInterface* visitor)
{
  _visitor = visitor;
  _deletevisitor = false;
}

const perulangan::IterativeSolverVisitorInterface* IterativeSolverWithVisitor::getVisitor() const
{
  assert(_visitor);
  return _visitor;
}

perulangan::IterativeSolverVisitorInterface* IterativeSolverWithVisitor::getVisitor()
{
  assert(_visitor);
  return _visitor;
}

const IterationInfo* IterativeSolverWithVisitor::getIterationInfo() const
{
  return &_info;
}

IterationInfo* IterativeSolverWithVisitor::getIterationInfo()
{
  return &_info;
}

alat::GhostVector& IterativeSolverWithVisitor::getMemory(int i) const
{
  if( i >= _memory.size() )
  {
    _error_string("getMemory", "", "too small memory");
  }
  return _memory[i];
}

/*--------------------------------------------------------------------------*/
void IterativeSolverWithVisitor::memory()
{
  _memory.set_size( getNVectors() );
  std::string type = getVisitor()->getVectorType();
  // int level = getVisitor()->getVectorLevel();
  // std::cerr << "IterativeSolverWithVisitor::memory() in " << getClassName() << " visitor = "<<getVisitor()->getClassName() << " " << getNVectors()  << " of type " << type << "\n";
  for(int i = 0; i < _memory.size(); i++)
  {
    std::stringstream ss;
    ss<<getClassName()<<"_"<<_id << "_" <<getVisitor()->getClassName()<<"_memory_"<<i;
    _memory[i] = alat::GhostVector( ss.str(), type);
    getVisitor()->newVector(&_memory[i]);
  }
}

/*--------------------------------------------------------------------------*/
// void IterativeSolverWithVisitor::basicInit(const alat::ParameterFile* parameterfile, std::string blockname)
// {
//   assert(_basicinitcalled == false);
//   // std::cerr << "IterativeSolver::basicInit() " << getClassName() << " blockname " << blockname <<"\n";
//   std::string paramblockname = blockname;
//   if(blockname == "none")
//   {
//     paramblockname = getClassName();
//   }
//   getVisitor()->basicInit(parameterfile, paramblockname);
//   getIterationInfo()->basicInit(parameterfile, paramblockname);
//   memory();
//   _basicinitcalled = true;
// }

#include  "Alat/agencylinearsolver.hpp"
#include  "Alat/tokenize.hpp"
#include  <cassert>

using namespace alat;

/*--------------------------------------------------------------------------*/
AgencyLinearSolver::~AgencyLinearSolver() {}
AgencyLinearSolver::AgencyLinearSolver() : alat::Map<alat::GhostLinearSolver, std::shared_ptr<perulangan::LinearSolverInterface> >() {}
AgencyLinearSolver::AgencyLinearSolver( const AgencyLinearSolver& vectoragency) : alat::Map<alat::GhostLinearSolver, std::shared_ptr<perulangan::LinearSolverInterface> >(vectoragency)
{
  assert(0);
}
AgencyLinearSolver& AgencyLinearSolver::operator=( const AgencyLinearSolver& vectoragency)
{
  alat::Map<alat::GhostLinearSolver, std::shared_ptr<perulangan::LinearSolverInterface> >::operator=(vectoragency);
  assert(0);
  return *this;
}
std::string AgencyLinearSolver::getClassName() const
{
  return "AgencyLinearSolver";
}
std::ostream& AgencyLinearSolver::printLoopInformation(std::ostream& os) const
{
  os << "AgencyLinearSolver: size=" << size() << "\n";
  for(const_iterator p=begin();p!=end();p++)
  {
    os << p->first << " --> " << p->second->getClassName() << "\n";
  }
  return os;
}

/*--------------------------------------------------------------------------*/
alat::StringIntMap AgencyLinearSolver::statistics() const
{
  alat::StringIntMap numberofdecriptions;
  for(const_iterator p = begin(); p != end(); p++)
  {
    std::string name = p->first.getClassName();
    alat::StringVector bouts = alat::Tokenize(name, "_");
    // name = "";
    // int last = 1;
    // if(bouts.size()>2) last = bouts.size()-2;
    // for(int i=0;i<last;i++)
    // {
    //   name += bouts[i]+"_";
    // }
    std::string typeandname = p->first.getType() + "_" + bouts[0];
    if( numberofdecriptions.find(typeandname) == numberofdecriptions.end() )
    {
      numberofdecriptions[typeandname] = 1;
    }
    else
    {
      numberofdecriptions[typeandname]++;
    }
  }
  return numberofdecriptions;
}

/*--------------------------------------------------------------------------*/
void AgencyLinearSolver::enrol(const alat::GhostLinearSolver& ghost)
{
  iterator p = find(ghost);
  if( p == end() )
  {
    ( *this )[ghost] = nullptr;//std::make_shared<perulangan::LinearSolverInterface>();
    // insert(std::make_pair(mg,static_cast<VectorInterface*>(NULL)));
  }
  else
  {
    std::cerr << "(register) already registered: " << ghost << std::endl;
    std::cerr << "(register) registered are:\n" << *this << std::endl;
    assert(0);
    exit(1);
  }
}

/*--------------------------------------------------------------------------*/
std::shared_ptr<const perulangan::LinearSolverInterface> AgencyLinearSolver::operator()(const alat::GhostLinearSolver& ghost) const
{
  const_iterator p = find(ghost);
  if( p == end() )
  {
    std::cerr << ": AgencyLinearSolver::operator(): ERROR"<<std::endl;
    std::cerr << ": alat::GhostLinearSolver '"<< ghost <<"' not found in list of: "<<std::endl;
    std::cerr << " "<< *this << std::endl;
    assert(0);
    exit(1);
  }
  std::shared_ptr<const perulangan::LinearSolverInterface>  vp = p->second;
  if(vp == nullptr)
  {
    std::cerr <<  "AgencyLinearSolver  std::shared_ptr<perulangan::LinearSolverInterface> is nullptr\t" << p->first;
    std::cerr << "\n" << *this << std::endl;
    assert(0);
    exit(1);
  }
  return vp;
}
std::shared_ptr<perulangan::LinearSolverInterface> AgencyLinearSolver::operator()(const alat::GhostLinearSolver& ghost)
{
  const_iterator p = find(ghost);
  if( p == end() )
  {
    std::cerr << ": AgencyLinearSolver::operator(): ERROR"<<std::endl;
    std::cerr << ": alat::GhostLinearSolver '"<< ghost <<"' not found in list of: "<<std::endl;
    std::cerr << " "<< *this << std::endl;
    assert(0);
    exit(1);
  }
  std::shared_ptr<perulangan::LinearSolverInterface>  vp = p->second;
  if(vp == nullptr)
  {
    std::cerr <<  "AgencyLinearSolver  std::shared_ptr<perulangan::LinearSolverInterface> is nullptr\t" << p->first;
    std::cerr << "\n" << *this << std::endl;
    assert(0);
    exit(1);
  }
  return vp;
}

/*--------------------------------------------------------------------------*/
std::ostream& alat::operator<<(std::ostream& os, const AgencyLinearSolver& vectoragency)
{
  int i = 0, n = vectoragency.size();
  os << "AgencyLinearSolver: size=" << n << ", ";
  for(alat::AgencyLinearSolver::const_iterator p = vectoragency.begin(); p != vectoragency.end(); p++, i++)
  {
    os << "GhostLinearSolver("<<i<<")=('"<< p->first << "',"<< p->second->getClassName() <<")\n";
    if( i < n-1 )
    {
      os << ", ";
    }
    else
    {
      os << ". ";
    }
  }
  return os;
}

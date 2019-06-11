#include  "Alat/agencymatrix.hpp"
#include  "Alat/tokenize.hpp"
#include  <cassert>

using namespace alat;

/*--------------------------------------------------------------------------*/
AgencyMatrix::~AgencyMatrix() {}
AgencyMatrix::AgencyMatrix() : alat::Map<alat::GhostMatrix, std::shared_ptr<alat::MatrixAllVariables> >() {}
AgencyMatrix::AgencyMatrix( const AgencyMatrix& matrixagency) : alat::Map<alat::GhostMatrix, std::shared_ptr<alat::MatrixAllVariables> >(matrixagency)
{
  assert(0);
}
AgencyMatrix& AgencyMatrix::operator=( const AgencyMatrix& matrixagency)
{
  alat::Map<alat::GhostMatrix, std::shared_ptr<alat::MatrixAllVariables> >::operator=(matrixagency);
  assert(0);
  return *this;
}
std::string AgencyMatrix::getClassName() const
{
  return "AgencyMatrix";
}
std::ostream& AgencyMatrix::printLoopInformation(std::ostream& os) const
{
  os << "AgencyMatrix: size=" << size() << "\n";
  for(const_iterator p=begin();p!=end();p++)
  {
    os << p->first << " --> " << p->second->getClassName() << "\n";
  }
  return os;
}

/*--------------------------------------------------------------------------*/
alat::StringIntMap AgencyMatrix::statistics() const
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
void AgencyMatrix::enrol(const alat::GhostMatrix& ghost)
{
  iterator p = find(ghost);
  if( p == end() )
  {
    ( *this )[ghost] = std::make_shared<alat::MatrixAllVariables>();
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
std::shared_ptr<const alat::MatrixAllVariables> AgencyMatrix::operator()(const alat::GhostMatrix& ghost) const
{
  const_iterator p = find(ghost);
  if( p == end() )
  {
    std::cerr << ": AgencyMatrix::operator(): ERROR"<<std::endl;
    std::cerr << ": alat::GhostMatrix '"<< ghost <<"' not found in list of: "<<std::endl;
    std::cerr << " "<< *this << std::endl;
    assert(0);
    exit(1);
  }
  std::shared_ptr<const alat::MatrixAllVariables>  vp = p->second;
  if(vp == nullptr)
  {
    std::cerr <<  "AgencyMatrix  std::shared_ptr<alat::MatrixAllVariables> is nullptr\t" << p->first;
    std::cerr << "\n" << *this << std::endl;
    assert(0);
    exit(1);
  }
  return vp;
}

/*--------------------------------------------------------------------------*/
std::shared_ptr<alat::MatrixAllVariables> AgencyMatrix::operator()(const alat::GhostMatrix& ghost)
{
  const_iterator p = find(ghost);
  if( p == end() )
  {
    std::cerr << ": AgencyMatrix::operator(): ERROR"<<std::endl;
    std::cerr << ": alat::GhostMatrix '"<< ghost <<"' not found in list of: "<<std::endl;
    std::cerr << " "<< *this << std::endl;
    assert(0);
    exit(1);
  }
  std::shared_ptr<alat::MatrixAllVariables>  vp = p->second;
  if(vp == nullptr)
  {
    std::cerr <<  "AgencyMatrix  std::shared_ptr<alat::MatrixAllVariables> is nullptr\t" << p->first;
    std::cerr << "\n" << *this << std::endl;
    assert(0);
    exit(1);
  }
  return vp;
}

/*--------------------------------------------------------------------------*/
std::ostream& alat::operator<<(std::ostream& os, const AgencyMatrix& matrixagency)
{
  int i = 0, n = matrixagency.size();
  os << "AgencyMatrix: size=" << n << ", ";
  for(alat::AgencyMatrix::const_iterator p = matrixagency.begin(); p != matrixagency.end(); p++, i++)
  {
    os << "GhostVector("<<i<<")=('"<< p->first << "',"<< p->second->getClassName() <<")\n";
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

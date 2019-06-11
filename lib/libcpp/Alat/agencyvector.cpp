#include  "Alat/agencyvector.hpp"
#include  "Alat/tokenize.hpp"
#include  <cassert>

using namespace alat;

/*--------------------------------------------------------------------------*/
AgencyVector::~AgencyVector() {}
AgencyVector::AgencyVector() : alat::Map<alat::GhostVector, std::shared_ptr<alat::VectorAllVariables> >() {}
AgencyVector::AgencyVector( const AgencyVector& vectoragency) : alat::Map<alat::GhostVector, std::shared_ptr<alat::VectorAllVariables> >(vectoragency)
{
  assert(0);
}
AgencyVector& AgencyVector::operator=( const AgencyVector& vectoragency)
{
  alat::Map<alat::GhostVector, std::shared_ptr<alat::VectorAllVariables> >::operator=(vectoragency);
  assert(0);
  return *this;
}
std::string AgencyVector::getClassName() const
{
  return "AgencyVector";
}
std::ostream& AgencyVector::printLoopInformation(std::ostream& os) const
{
  os << "AgencyVector: size=" << size() << "\n";
  for(const_iterator p=begin();p!=end();p++)
  {
    os << p->first << " --> " << p->second->getClassName() << "\n";
  }
  return os;
}

/*--------------------------------------------------------------------------*/
alat::StringIntMap AgencyVector::statistics() const
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
void AgencyVector::enrol(const alat::GhostVector& ghost)
{
  iterator p = find(ghost);
  if( p == end() )
  {
    ( *this )[ghost] = std::make_shared<alat::VectorAllVariables>();
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
std::shared_ptr<const alat::VectorAllVariables> AgencyVector::operator()(const alat::GhostVector& ghost) const
{
  const_iterator p = find(ghost);
  if( p == end() )
  {
    std::cerr << ": AgencyVector::operator(): ERROR"<<std::endl;
    std::cerr << ": alat::GhostVector '"<< ghost <<"' not found in list of: "<<std::endl;
    std::cerr << " "<< *this << std::endl;
    assert(0);
    exit(1);
  }
  std::shared_ptr<const alat::VectorAllVariables>  vp = p->second;
  if(vp == nullptr)
  {
    std::cerr <<  "AgencyVector  std::shared_ptr<alat::VectorAllVariables> is nullptr\t" << p->first;
    std::cerr << "\n" << *this << std::endl;
    assert(0);
    exit(1);
  }
  return vp;
}
std::shared_ptr<alat::VectorAllVariables> AgencyVector::operator()(const alat::GhostVector& ghost)
{
  const_iterator p = find(ghost);
  if( p == end() )
  {
    std::cerr << ": AgencyVector::operator(): ERROR"<<std::endl;
    std::cerr << ": alat::GhostVector '"<< ghost <<"' not found in list of: "<<std::endl;
    std::cerr << " "<< *this << std::endl;
    assert(0);
    exit(1);
  }
  std::shared_ptr<alat::VectorAllVariables>  vp = p->second;
  if(vp == nullptr)
  {
    std::cerr <<  "AgencyVector  std::shared_ptr<alat::VectorAllVariables> is nullptr\t" << p->first;
    std::cerr << "\n" << *this << std::endl;
    assert(0);
    exit(1);
  }
  return vp;
}

/*--------------------------------------------------------------------------*/
std::ostream& alat::operator<<(std::ostream& os, const AgencyVector& vectoragency)
{
  int i = 0, n = vectoragency.size();
  os << "AgencyVector: size=" << n << ", ";
  for(alat::AgencyVector::const_iterator p = vectoragency.begin(); p != vectoragency.end(); p++, i++)
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

#include  "Alat/strings.hpp"
#include  "Alat/stringvector.hpp"
#include  "Alat/tokenize.hpp"
#include  <iterator>

using namespace alat;

/*--------------------------------------------------------------------------*/
StringVector::~StringVector() {}
StringVector::StringVector() : alat::Vector<std::string>() {}
StringVector::StringVector(const StringVector& v) : alat::Vector<std::string>(v)  {}
StringVector::StringVector(const StringSet& v) : alat::Vector<std::string>()
{
  set_size( v.size() );
  std::copy( v.begin(), v.end(), begin() );
}

StringVector::StringVector(int n) : alat::Vector<std::string>(n, "") {}
StringVector::StringVector(std::string valuechain, const std::string& sep) : alat::Vector<std::string>()
{
  reInitFromString(valuechain, sep);
}

StringVector::StringVector(std::istream_iterator<std::string> begin, std::istream_iterator<std::string> end) : alat::Vector<std::string>(begin, end) {}

StringVector& StringVector::operator=(std::string n)
{
  alat::StringVector::equal(n);
  return *this;
}

StringVector& StringVector::operator=(const StringVector& v)
{
  alat::StringVector::equal(v);
  return *this;
}

/*--------------------------------------------------------------------------*/
void StringVector::reInitFromString(std::string valuechain, const std::string& sep)
{
  assert(0);
  alat::StringVector bouts = alat::Tokenize(valuechain, sep);
  int n = bouts.size();
  Vector<std::string>::set_size(n);
  for(int i = 0; i < n; i++)
  {
    ( *this )[i] =  bouts[i].c_str();
  }
}

/*--------------------------------------------------------------------------*/
bool StringVector::contains(const StringSet& v) const
{
  if( size() != v.size() )
  {
    return false;
  }
  for(const_iterator p = begin(); p != end(); p++)
  {
    if( v.find(*p) == v.end() )
    {
      return false;
    }
  }
  return true;
}

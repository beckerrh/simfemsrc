#include  "Alat/strings.hpp"
#include  "Alat/tokenize.hpp"
#include  <iterator>

using namespace alat;

/*--------------------------------------------------------------------------*/
StringSet::~StringSet() {}
StringSet::StringSet() : alat::Set<std::string>() {}
StringSet::StringSet(const StringSet& v) : alat::Set<std::string>(v)  {}
StringSet::StringSet(const StringVector& v)
{
  clear();
  for(int i=0;i<v.size();i++)
  {
    insert(v[i]);
  }
}
StringSet::StringSet(std::string valuechain, const std::string& sep) : alat::Set<std::string>()
{
  reInitFromString(valuechain, sep);
}
StringSet::StringSet(const char* valuechain, const std::string& sep) : alat::Set<std::string>()
{
  reInitFromString(std::string(valuechain), sep);
}

StringSet& StringSet::operator=(const StringSet& v)
{
  alat::Set<std::string>::operator=(v);
  return *this;
}
void StringSet::reInitFromString(std::string valuechain, const std::string& sep)
{
  alat::StringVector bouts = alat::Tokenize(valuechain, sep);
  for(int i = 0; i < bouts.size(); i++)
  {
    insert(bouts[i].c_str());
  }
}

/*--------------------------------------------------------------------------*/
StringList::~StringList() {}
StringList::StringList() : alat::List<std::string>() {}
StringList::StringList(const StringList& v) : alat::List<std::string>(v)  {}
StringList::StringList(const StringVector& v)
{
  clear();
  for(int i=0;i<v.size();i++)
  {
    push_back(v[i]);
  }
}
StringList::StringList(std::string valuechain, const std::string& sep) : alat::List<std::string>()
{
  reInitFromString(valuechain, sep);
}
StringList::StringList(const char* valuechain, const std::string& sep) : alat::List<std::string>()
{
  reInitFromString(std::string(valuechain), sep);
}

StringList& StringList::operator=(const StringList& v)
{
  alat::List<std::string>::operator=(v);
  return *this;
}
void StringList::reInitFromString(std::string valuechain, const std::string& sep)
{
  alat::StringVector bouts = alat::Tokenize(valuechain, sep);
  for(int i = 0; i < bouts.size(); i++)
  {
    push_back(bouts[i].c_str());
  }
}

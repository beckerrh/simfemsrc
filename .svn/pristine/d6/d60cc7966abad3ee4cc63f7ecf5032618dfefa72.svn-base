#ifndef __Alat_Strings_h
#define __Alat_Strings_h

#include  "set.hpp"
#include  "list.hpp"

/*---------------------------------------------------------*/

namespace alat
{
  class StringVector;

  class StringSet : public alat::Set<std::string>
  {
  public:
    ~StringSet();
    StringSet();
    StringSet(const StringSet& v);
    StringSet(const StringVector& v);
    StringSet(const char* valuechain, const std::string& sep=",");
    StringSet(std::string valuechain, const std::string& sep=",");
    void reInitFromString(std::string valuechain, const std::string& sep=",");
    StringSet& operator=(const StringSet& v);
  };

  class StringList : public alat::List<std::string>
  {
  public:
    ~StringList();
    StringList();
    StringList(const StringList& v);
    StringList(const StringVector& v);
    StringList(const char* valuechain, const std::string& sep=",");
    StringList(std::string valuechain, const std::string& sep=",");
    void reInitFromString(std::string valuechain, const std::string& sep=",");
    StringList& operator=(const StringList& v);
  };
}

/*---------------------------------------------------------*/

#endif

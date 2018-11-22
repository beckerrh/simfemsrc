#ifndef __alat_Ghost_h
#define __alat_Ghost_h

#include  <string>

/*--------------------------------------------------------------------------*/
namespace alat
{
  class Ghost
  {
private:
    std::string _name, _type;

public:
    virtual ~Ghost();
    Ghost();
    Ghost( const Ghost& ghost);
    Ghost(const std::string name);
    Ghost(const std::string name, const std::string type);
    Ghost& operator=( const Ghost& ghost);

    void setName(const std::string& name);
    const std::string& getName() const;
    void setType(const std::string& type);
    const std::string& getType() const;
    virtual std::string getClassName() const;
    bool operator==(const Ghost& v) const;
    bool operator<(const Ghost& v) const;
  };
  std::ostream& operator<<(std::ostream& os, const Ghost& g);
}

/*--------------------------------------------------------------------------*/

#endif

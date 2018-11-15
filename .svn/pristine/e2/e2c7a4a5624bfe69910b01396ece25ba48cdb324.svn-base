#ifndef __alat_GhostVector_h
#define __alat_GhostVector_h

#include  "Alat/ghost.hpp"

/*--------------------------------------------------------------------------*/

namespace alat
{
  class GhostVector : public Ghost
  {
private:

public:
    ~GhostVector();
    GhostVector();
    GhostVector( const GhostVector& ghost);
    GhostVector(const std::string name);
    GhostVector(const std::string name, const std::string type);
    std::string getClassName() const;
    GhostVector& operator=( const GhostVector& ghost);
  };
  std::ostream& operator<<(std::ostream& os, const GhostVector& g);
}

/*--------------------------------------------------------------------------*/

#endif

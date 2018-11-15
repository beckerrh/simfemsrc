#ifndef __alat_GhostMatrix_h
#define __alat_GhostMatrix_h

#include  "Alat/ghost.hpp"

/*--------------------------------------------------------------------------*/

namespace alat
{
  class GhostMatrix : public Ghost
  {
private:

public:
    ~GhostMatrix();
    GhostMatrix();
    GhostMatrix( const GhostMatrix& ghost);
    GhostMatrix(const std::string name);
    GhostMatrix(const std::string name, const std::string description);
    GhostMatrix& operator=( const GhostMatrix& ghost);
    std::string getClassName() const;
  };
  std::ostream& operator<<(std::ostream& os, const GhostMatrix& g);
}

/*--------------------------------------------------------------------------*/

#endif

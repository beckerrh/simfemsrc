#ifndef __Perulangan_RichardsonOptimal_h
#define __Perulangan_RichardsonOptimal_h

#include  "richardsonrb.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class RichardsonOptimal : public RichardsonRB
  {
private:
    void _computeSmallSystem(int index, int nmemory) const;
    void _addvector(alat::GhostVector& u, int nmemory) const;
    void _matrixVectorProduct(int index) const;

public:
    ~RichardsonOptimal();
    RichardsonOptimal(const std::string& type, int nvectors, const std::string& solutiontype);
    RichardsonOptimal( const RichardsonOptimal& richardsonoptimal);
    RichardsonOptimal& operator=( const RichardsonOptimal& richardsonoptimal);
    std::string getClassName() const;

    int getNVectors() const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

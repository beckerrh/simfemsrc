#ifndef __Perulangan_RichardsonSystem_h
#define __Perulangan_RichardsonSystem_h

#include  "richardsonrb.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class RichardsonSystem : public RichardsonRB
  {
private:
    mutable arma::vec _scalarproduct;
    mutable alat::Vector<alat::GhostVector*> _vectors;

    void _computeSmallSystem(int index, int nmemory) const;
    void _addvector(alat::GhostVector& u, int nmemory) const;
    void _matrixVectorProduct(int index) const;

public:
    ~RichardsonSystem();
    RichardsonSystem(const std::string& type, int nvectors, const std::string& solutiontype);
    RichardsonSystem( const RichardsonSystem& richardsonoptimal);
    RichardsonSystem& operator=( const RichardsonSystem& richardsonoptimal);
    std::string getClassName() const;

    void basicInit(const alat::ParameterFile* parameterfile, std::string blockname);
    int getNVectors() const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

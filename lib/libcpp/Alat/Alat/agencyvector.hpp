#ifndef __solvers_AgencyVector_h
#define __solvers_AgencyVector_h

#include  "Alat/ghostvector.hpp"
#include  "Alat/map.hpp"
#include  "Alat/vectorallvariables.hpp"

/*--------------------------------------------------------------------------*/

namespace alat
{
  class AgencyVector : public alat::Map<alat::GhostVector, std::shared_ptr<alat::VectorAllVariables> >
  {
public:
    typedef alat::Map<alat::GhostVector, std::shared_ptr<alat::VectorAllVariables> >::const_iterator const_iterator;
    typedef alat::Map<alat::GhostVector, std::shared_ptr<alat::VectorAllVariables> >::iterator iterator;

public:
    ~AgencyVector();
    AgencyVector();
    AgencyVector( const AgencyVector& vectoragency);
    AgencyVector& operator=( const AgencyVector& vectoragency);
    std::string getClassName() const;
    std::ostream& printLoopInformation(std::ostream& os) const;
    void enrol(const alat::GhostVector& ghost);
    std::shared_ptr<const alat::VectorAllVariables> operator()(const alat::GhostVector& ghost) const;
    std::shared_ptr<alat::VectorAllVariables> operator()(const alat::GhostVector& ghost);
    alat::StringIntMap statistics() const;
  };
  std::ostream& operator<<(std::ostream& os, const AgencyVector& vectoragency);
}

/*--------------------------------------------------------------------------*/

#endif

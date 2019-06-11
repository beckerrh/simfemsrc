#ifndef __solvers_AgencyMatrix_h
#define __solvers_AgencyMatrix_h

#include  "Alat/ghostmatrix.hpp"
#include  "Alat/map.hpp"
#include  "Alat/matrixallvariables.hpp"

/*--------------------------------------------------------------------------*/

namespace alat
{
  class AgencyMatrix : public alat::Map<alat::GhostMatrix, std::shared_ptr<alat::MatrixAllVariables> >
  {
public:
    typedef alat::Map<alat::GhostMatrix, std::shared_ptr<alat::MatrixAllVariables> >::const_iterator const_iterator;
    typedef alat::Map<alat::GhostMatrix, std::shared_ptr<alat::MatrixAllVariables> >::iterator iterator;

public:
    ~AgencyMatrix();
    AgencyMatrix();
    AgencyMatrix( const AgencyMatrix& matrixagency);
    AgencyMatrix& operator=( const AgencyMatrix& matrixagency);
    std::string getClassName() const;
    std::ostream& printLoopInformation(std::ostream& os) const;
    void enrol(const alat::GhostMatrix& ghost);
    std::shared_ptr<const alat::MatrixAllVariables> operator()(const alat::GhostMatrix& ghost) const;
    std::shared_ptr<alat::MatrixAllVariables> operator()(const alat::GhostMatrix& ghost);
    alat::StringIntMap statistics() const;
  };
  std::ostream& operator<<(std::ostream& os, const AgencyMatrix& matrixagency);
}

/*--------------------------------------------------------------------------*/

#endif

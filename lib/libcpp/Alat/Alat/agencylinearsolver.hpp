#ifndef __solvers_AgencyLinearSolver_h
#define __solvers_AgencyLinearSolver_h

#include  "Alat/ghostlinearsolver.hpp"
#include  "Alat/map.hpp"
#include  "Alat/vectorallvariables.hpp"
#include  "Perulangan/linearsolverinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class AgencyLinearSolver : public alat::Map<alat::GhostLinearSolver, std::shared_ptr<perulangan::LinearSolverInterface> >
  {
public:
    typedef alat::Map<alat::GhostLinearSolver, std::shared_ptr<perulangan::LinearSolverInterface> >::const_iterator const_iterator;
    typedef alat::Map<alat::GhostLinearSolver, std::shared_ptr<perulangan::LinearSolverInterface> >::iterator iterator;

public:
    ~AgencyLinearSolver();
    AgencyLinearSolver();
    AgencyLinearSolver( const AgencyLinearSolver& vectoragency);
    AgencyLinearSolver& operator=( const AgencyLinearSolver& vectoragency);
    std::string getClassName() const;
    std::ostream& printLoopInformation(std::ostream& os) const;
    void enrol(const alat::GhostLinearSolver& ghost);
    std::shared_ptr<const perulangan::LinearSolverInterface> operator()(const alat::GhostLinearSolver& ghost) const;
    std::shared_ptr<perulangan::LinearSolverInterface> operator()(const alat::GhostLinearSolver& ghost);
    alat::StringIntMap statistics() const;
  };
  std::ostream& operator<<(std::ostream& os, const AgencyLinearSolver& vectoragency);
}

/*--------------------------------------------------------------------------*/

#endif

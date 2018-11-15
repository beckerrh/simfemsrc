#ifndef __alat_GhostLinearSolver_h
#define __alat_GhostLinearSolver_h

#include  "Alat/ghostmatrix.hpp"

/*--------------------------------------------------------------------------*/

namespace alat
{
  class GhostLinearSolver : public Ghost
  {
private:
   const alat::GhostMatrix _matrix;

public:
    ~GhostLinearSolver();
    GhostLinearSolver( const GhostLinearSolver& ghost);
    GhostLinearSolver(const std::string name, const std::string description, const alat::GhostMatrix matrix);
    GhostLinearSolver& operator=( const GhostLinearSolver& ghost);
    std::string getClassName() const;
    const alat::GhostMatrix& getMatrix() const;
  };
  std::ostream& operator<<(std::ostream& os, const GhostLinearSolver& g);
}

/*--------------------------------------------------------------------------*/

#endif

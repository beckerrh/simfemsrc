#include  "Perulangan/enums.hpp"
#include  <cassert>
#include  <iostream>

namespace perulanganEnums
{
  std::string iterationStatusToString(iterationstatus i)
  {
    if(i == IterationStatusNone)
    {
      return "IterationStatusNone";
    }
    else if(i == IterationStatusConverged)
    {
      return "IterationStatusConverged";
    }
    else if(i == IterationStatusDiverged)
    {
      return "IterationStatusDiverged";
    }
    else if(i == IterationStatusMaxIter)
    {
      return "IterationStatusMaxIter";
    }
    else if(i == IterationStatusRunning)
    {
      return "IterationStatusRunning";
    }
    else if(i == IterationStatusProblem)
    {
      return "IterationStatusProblem";
    }
    std::cerr<<"****iterationStatusToString: iterationstatus not defined "<<i<<'\n';
    assert(0);
  }

  std::string residualStatusToString(residualstatus r)
  {
    if(r == ResidualStatusNone)
    {
      return "ResidualStatusNone";
    }
    else if(r == ResidualStatusOk)
    {
      return "ResidualStatusOk";
    }
    else if(r == ResidualStatusNotOk)
    {
      return "ResidualStatusNotOk";
    }
    std::cerr<<"****residualStatusToString: residualstatus not defined "<<r<<'\n';
    assert(0);
  }

  std::string matrixStatusToString(matrixstatus r)
  {
    if(r == MatrixStatusNone)
    {
      return "MatrixStatusNone";
    }
    else if(r == MatrixStatusOk)
    {
      return "MatrixStatusOk";
    }
    else if(r == MatrixStatusNotOk)
    {
      return "MatrixStatusNotOk";
    }
    std::cerr<<"****matrixStatusToString: matrixstatus not defined "<<r<<'\n';
    assert(0);
  }

  std::string newtonStatusToString(newtonstatus n)
  {
    if(n == NewtonStatusNone)
    {
      return "NewtonStatusNone";
    }
    else if(n == NewtonStatusResidualStatusNotOk)
    {
      return "NewtonStatusResidualStatusNotOk";
    }
    else if(n == NewtonStatusConverged)
    {
      return "NewtonStatusConverged";
    }
    else if(n == NewtonStatusRunning)
    {
      return "NewtonStatusRunning";
    }
    else if(n == NewtonStatusLinearNotOk)
    {
      return "NewtonStatusLinearNotOk";
    }
    else if(n == NewtonStatusBadReduction)
    {
      return "NewtonStatusBadReduction";
    }
    else if(n == NewtonStatusTooManyIterations)
    {
      return "NewtonStatusTooManyIterations";
    }
    else if(n == NewtonStatusDiverged)
    {
      return "NewtonStatusDiverged";
    }
    else if(n == NewtonStatusMaxLineSearchAttained)
    {
      return "NewtonStatusMaxLineSearchAttained";
    }
    std::cerr<<"****newtonStatusToString: newtonstatus not defined "<<n<<'\n';
    assert(0);
  }
}

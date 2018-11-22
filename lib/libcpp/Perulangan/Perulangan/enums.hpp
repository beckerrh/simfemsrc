#ifndef  __Perulangan_enum_h
#define  __Perulangan_enum_h

#include  <string>

/*---------------------------------------------*/
namespace perulanganEnums
{
  enum iterationstatus {IterationStatusNone, IterationStatusConverged, IterationStatusDiverged, IterationStatusMaxIter, IterationStatusRunning, IterationStatusProblem};
  enum matrixstatus {MatrixStatusNone, MatrixStatusOk, MatrixStatusNotOk};
  enum residualstatus {ResidualStatusNone, ResidualStatusOk, ResidualStatusNotOk};
  enum newtonstatus {NewtonStatusNone,NewtonStatusResidualStatusNotOk, NewtonStatusConverged, NewtonStatusDiverged, NewtonStatusRunning, NewtonStatusLinearNotOk,  NewtonStatusBadReduction, NewtonStatusTooManyIterations, NewtonStatusMaxLineSearchAttained};

  std::string iterationStatusToString(iterationstatus i);
  std::string matrixStatusToString(matrixstatus r);
  std::string residualStatusToString(residualstatus r);
  std::string newtonStatusToString(newtonstatus n);
}

#endif

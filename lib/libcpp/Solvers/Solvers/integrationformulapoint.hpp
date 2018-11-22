#ifndef __Solvers_IntegrationFormulaPoint_h
#define __Solvers_IntegrationFormulaPoint_h

#include  "Solvers/integrationformulainterface.hpp"

/*---------------------------------------------------------*/
namespace solvers
{
  /*------------------------------------------------------------*/
  class IntegrationFormulaPoint : public IntegrationFormulaInterface
  {
  public:
    std::string getClassName() const{return "IntegrationFormulaPoint";}
    IntegrationFormulaPoint();
  };
  /*---------------------------------------------------------*/
}
#endif

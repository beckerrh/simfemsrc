#include  "Solvers/integrationformulapoint.hpp"

using namespace solvers;
using namespace std;

/*------------------------------------------------------------*/
IntegrationFormulaPoint::IntegrationFormulaPoint() : IntegrationFormulaInterface(1)
{
  _w[0] = 1.0;
  _c[0].x() = 0.;
}

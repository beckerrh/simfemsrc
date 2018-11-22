#include  "Solvers/integrationformulatetrahedral.hpp"

using namespace solvers;

/*--------------------------------------------------------------------------*/
TetrahedralFormula1::TetrahedralFormula1() : IntegrationFormulaInterface(1)
{
  _w[0] = 1.0;
  _c[0].x() = 0.25;
  _c[0].y() = 0.25;
  _c[0].z() = 0.25;
}
std::string TetrahedralFormula1::getClassName() const {return "TetrahedralFormula1";}

/*--------------------------------------------------------------------------*/
TetrahedronTrapez::TetrahedronTrapez() : IntegrationFormulaInterface(4)
{
  _w[0] = 0.25;
  _w[1] = _w[0];
  _w[2] = _w[0];
  _w[3] = _w[0];

  _c[0].x() = 0.0;
  _c[0].y() = 0.0;
  _c[0].z() = 0.0;
  _c[1].x() = 1.0;
  _c[1].y() = 0.0;
  _c[1].z() = 0.0;
  _c[2].x() = 0.0;
  _c[2].y() = 1.0;
  _c[2].z() = 0.0;
  _c[3].x() = 0.0;
  _c[3].y() = 0.0;
  _c[3].z() = 1.0;
}
std::string TetrahedronTrapez::getClassName() const {return "TetrahedronTrapez";}

/*--------------------------------------------------------------------------*/
TetrahedralFormula4::TetrahedralFormula4() : IntegrationFormulaInterface(4)
{
  _w[0] = 0.25;
  _w[1] = _w[0];
  _w[2] = _w[0];
  _w[3] = _w[0];

  const double a = 0.585410196624969;
  const double b = 0.138196601125011;

  _c[0].x() = a;
  _c[0].y() = b;
  _c[0].z() = b;
  _c[1].x() = b;
  _c[1].y() = a;
  _c[1].z() = b;
  _c[2].x() = b;
  _c[2].y() = b;
  _c[2].z() = a;
  _c[3].x() = b;
  _c[3].y() = b;
  _c[3].z() = b;
}
std::string TetrahedralFormula4::getClassName() const {return "TetrahedralFormula4";}

/*--------------------------------------------------------------------------*/
TetrahedralFormula5::TetrahedralFormula5() : IntegrationFormulaInterface(5)
{
  // degree 3
  _w[0] = 0.444444444444444;
  _w[1] = 0.138888886888889;
  _w[2] = _w[1];
  _w[3] = _w[1];
  _w[4] = _w[1];
  const double a = 0.724086765841831;
  const double b = 0.091971078062723;
  _c[0].x() = 0.25;
  _c[0].y() = 0.25;
  _c[0].z() = 0.25;

  _c[1].x() = a;
  _c[1].y() = b;
  _c[1].z() = b;

  _c[2].x() = b;
  _c[2].y() = a;
  _c[2].z() = b;

  _c[3].x() = b;
  _c[3].y() = b;
  _c[3].z() = a;

  _c[4].x() = b;
  _c[4].y() = b;
  _c[4].z() = b;
}
std::string TetrahedralFormula5::getClassName() const {return "TetrahedralFormula5";}

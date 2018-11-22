#include  <string>
#include  <cassert>
#include  <iostream>
#include  <boost/python.hpp>
#include  <boost/python/numpy.hpp>
#include  "Tools/armatopy.hpp"
#include  "solverwrapper.hpp"

namespace bp = boost::python;
namespace np = boost::python::numpy;

void (solvers::Solver::*sp1)(std::string, std::string)  = &solvers::Solver::setParameter;
void (solvers::Solver::*sp2)(std::string, int)  = &solvers::Solver::setParameter;
void (solvers::Solver::*sp3)(std::string, double)  = &solvers::Solver::setParameter;
void (solvers::Solver::*sp4)(std::string, bool)  = &solvers::Solver::setParameter;

/*--------------------------------------------------------------------------*/
void wrapSolver()
{
  bp::to_python_converter<solvers::ErrorsMap, simfem::ErrorsMapToDictConverter>();
  bp::to_python_converter<arma::vec, simfem::ArmavecToNumpyConverter>();
  bp::class_< solvers::Solver >("Solver")
  .def("loadMesh", &solvers::Solver::loadMesh)
  .def("setMesh", &solvers::Solver::setMesh)
  .def("init", &solvers::Solver::init)
  .def("writeXdmf", &solvers::Solver::writeXdmf)
  // .def("addVariablePlain", &solvers::Solver::addVariablePlain)
  .def("setParameter", sp1)
  .def("setParameter", sp2)
  .def("setParameter", sp3)
  .def("setParameter", sp4)
  .def("getInfo", &solvers::Solver::getInfo)
  .def("getMeshInfo", &solvers::Solver::getMeshInfo)
  .def(bp::self_ns::str(bp::self_ns::self));
}

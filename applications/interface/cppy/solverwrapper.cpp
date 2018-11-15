#include  <string>
#include  <cassert>
#include  <iostream>
#include  <boost/python.hpp>
#include  <boost/python/numpy.hpp>
#include  "solver.hpp"
#include  "Tools/armatopy.hpp"
#include  "MeshWrap/meshwrapper.hpp"
#include  "SolversWrap/solverwrapper.hpp"

namespace bp = boost::python;
namespace np = boost::python::numpy;

/*--------------------------------------------------------------------------*/
void wrapSolver()
{
  bp::class_<Solver, bp::bases<solvers::Solver> >("Solver")
  // overwritten functions
  .def("init", &Solver::init)
  .def("run", &Solver::run)
  .def(bp::self_ns::str(bp::self_ns::self));
}

/*---------------------------------------------------------------------------*/
BOOST_PYTHON_MODULE(simfeminterface)
{
  Py_Initialize();
  np::initialize();
  wrapSolver();
}

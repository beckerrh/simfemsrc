#include  <string>
#include  <cassert>
#include  <iostream>
#include  <boost/python.hpp>
#include  <boost/python/numpy.hpp>
#include  "MeshWrap/meshwrapper.hpp"
#include  "MeshWrap/timemeshwrapper.hpp"
#include  "SolversWrap/solverwrapper.hpp"

namespace bp = boost::python;
namespace np = boost::python::numpy;

/*---------------------------------------------------------------------------*/
BOOST_PYTHON_MODULE(simfempy)
{
  Py_Initialize();
  np::initialize();

  wrapMesh();
  wrapTimeMesh();
  wrapSolver();
}

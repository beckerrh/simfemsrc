#include  <string>
#include  <cassert>
#include  <iostream>
#include  <boost/python.hpp>
#include  <boost/python/numpy.hpp>
#include  "timemeshwrapper.hpp"

namespace bp = boost::python;
namespace np = boost::python::numpy;

/*---------------------------------------------------------------------------*/
void wrapTimeMesh()
{
  bp::class_< TimeMeshWrapper >("TimeMesh")
  .def(bp::init<int, double, double>())
  .def("setData", &TimeMeshWrapper::setData)
  .def(bp::self_ns::str(bp::self_ns::self));
}

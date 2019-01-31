#include  "Solvers/solver.hpp"
#include  <boost/python.hpp>
#include  <boost/python/numpy.hpp>
#include  "Alat/armadillo.hpp"

namespace bp = boost::python;
namespace np = boost::python::numpy;

/*---------------------------------------------------------------------------*/

namespace simfem
{
  struct ArmavecToNumpyConverter{
    static PyObject* convert(const alat::armavec& avec);
  };
  struct ErrorsMapToDictConverter{
    static PyObject* convert(const solvers::ErrorsMap& map);
  };

  alat::armavec numpy_to_armavec(const np::ndarray& y);
  boost::python::object armavec_to_numpy(const alat::armavec& v);
  arma::mat numpy_to_armamat(const np::ndarray& y);
  boost::python::object armamat_to_numpy(const arma::mat& v);
}

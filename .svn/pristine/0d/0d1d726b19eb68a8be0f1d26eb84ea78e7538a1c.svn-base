#include  "Tools/armatopy.hpp"
#include  <cassert>
#include  <iostream>

namespace bp = boost::python;
namespace np = boost::python::numpy;
using namespace simfem;

/*---------------------------------------------------------------------------*/
PyObject* ArmavecToNumpyConverter::convert(const arma::vec& avec)
{
  // std::cerr << "@@@@@@@@-------------- " << avec << "\n";
  const double* data = avec.memptr();
  // std::cerr << "nr nc =" << avec.n_rows << " " << avec.n_cols << "\n";
  bp::tuple shape = bp::make_tuple(avec.n_rows);
  // std::cerr << "arma_to_numpy() " << bp::extract<char const *>(bp::str(shape)) << "\n";
  bp::tuple stride = bp::make_tuple(sizeof(double));
  np::dtype dt = np::dtype::get_builtin<double>();
  // np::ndarray nparray = np::from_data(data, dt, shape, stride, bp::object());
  bp::object own;
  np::ndarray nparray = np::from_data(data, dt, shape, stride, own);
  // std::cerr << "array created" << std::endl<< bp::extract<char const *>(bp::str(nparray)) << std::endl;
  return bp::incref(nparray.copy().ptr());
  // return incref(make_tuple(pair.first, pair.second).ptr());
}
/*---------------------------------------------------------------------------*/
PyObject* ErrorsMapToDictConverter::convert(const solvers::ErrorsMap& map) {
  bp::dict pmap;
  for(solvers::ErrorsMap::const_iterator p = map.begin();p!=map.end();p++)
  {
    boost::python::object pyvec(boost::python::handle<>(simfem::ArmavecToNumpyConverter::convert(p->second)));
    pmap[p->first] = pyvec;
  }
  return bp::incref(pmap.ptr());
}

/*---------------------------------------------------------------------------*/
arma::vec simfem::numpy_to_armavec(const np::ndarray& y)
{
  // std::cerr << "shape="<<y.shape(0)<< " " << y.shape(1)<<"\n";
  return arma::vec(reinterpret_cast<double*>(y.get_data()),y.shape(0));
}
boost::python::object simfem::armavec_to_numpy(const arma::vec& v)
{
  // std::cerr << "arma_to_numpy()\n";
  // int data[] = {1,2,3,4,5};
  // v.print("hello");
  const double* data = v.memptr();
  // std::cerr << "nr nc =" << v.n_rows << " " << v.n_cols << "\n";
  bp::tuple shape = bp::make_tuple(v.n_rows);
  // std::cerr << "arma_to_numpy() " << bp::extract<char const *>(bp::str(shape)) << "\n";
  bp::tuple stride = bp::make_tuple(sizeof(double));
  np::dtype dt = np::dtype::get_builtin<double>();
  np::ndarray nparray = np::from_data(data, dt, shape, stride, bp::object());
  // std::cerr << "nparray=" << bp::extract<char const *>(bp::str(nparray)) << "\n";
  return nparray;
}
arma::mat simfem::numpy_to_armamat(const np::ndarray& y)
{
  // std::cerr << "shape="<<y.shape(0)<< " " << y.shape(1)<<"\n";
  // return arma::mat(reinterpret_cast<double*>(y.get_data()),y.shape(0),y.shape(1));
  return arma::mat(reinterpret_cast<double*>(y.get_data()),y.shape(0)/2,2);
}
boost::python::object simfem::armamat_to_numpy(const arma::mat& v)
{
  const double* data = v.memptr();
  // std::cerr << "nr nc =" << v.n_rows << " " << v.n_cols << "\n";
  // bp::tuple shape = bp::make_tuple(v.n_cols, v.n_rows);
  // bp::tuple shape = bp::make_tuple(v.n_rows, v.n_cols);
  bp::tuple shape = bp::make_tuple(v.n_rows*v.n_cols);
  // std::cerr << "armamat_to_numpy() " << bp::extract<char const *>(bp::str(shape)) << "\n";
  // bp::tuple stride = bp::make_tuple(sizeof(double),v.n_rows*sizeof(double));
  bp::tuple stride = bp::make_tuple(sizeof(double));
  np::dtype dt = np::dtype::get_builtin<double>();
  // std::cerr << "armamat_to_numpy() " << bp::extract<char const *>(bp::str(stride)) << "\n";
  np::ndarray nparray = np::from_data(data, dt, shape, stride, bp::object());
  // std::cerr << "nparray=" << bp::extract<char const *>(nparray) << "\n";
  return nparray;
}

#include  <string>
#include  <cassert>
#include  <iostream>
#include  <boost/python.hpp>
#include  <boost/python/numpy.hpp>
#include  "meshwrapper.hpp"
// #include  "Mesh/mesh.hpp"

/*---------------------------------------------------------------------------*/
template<class T> T * get_pointer( std::shared_ptr<T> const& p ) {return p.get();}
/*---------------------------------------------------------------------------*/

namespace bp = boost::python;
namespace np = boost::python::numpy;

// /*---------------------------------------------------------------------------*/
// void wrapMesh()
// {
//   bp::def("create", create);
// 	bp::class_< mesh::Mesh >("Mesh", bp::no_init)
// 	.def("addGeometryObject", &mesh::Mesh::addGeometryObjectByName)
// 	.def("readGmsh", &mesh::Mesh::readGmsh)
// 	.def("writeVtk", &mesh::Mesh::writeVtk)
// 	.def("writeBoundaryVtk", &mesh::Mesh::writeBoundaryVtk)
// 	.def("load", &mesh::Mesh::loadH5)
// 	.def("save", &mesh::Mesh::saveH5)
// 	.def(bp::self_ns::str(bp::self_ns::self));
//   bp::register_ptr_to_python< std::shared_ptr<mesh::Mesh> >();
//   bp::implicitly_convertible<std::shared_ptr<mesh::Mesh>,std::shared_ptr<mesh::MeshInterface> >();
// }
//
/*---------------------------------------------------------------------------*/
void wrapMesh()
{
  bp::def("create", create);
	bp::class_< MeshWrapper >("Mesh", bp::no_init)
	.def("addGeometryObject", &MeshWrapper::addGeometryObjectByName)
	.def("readGmsh", &MeshWrapper::readGmsh)
	.def("writeVtk", &MeshWrapper::writeVtk)
	.def("writeBoundaryVtk", &MeshWrapper::writeBoundaryVtk)
	.def("load", &MeshWrapper::loadH5)
	.def("save", &MeshWrapper::saveH5)
  .def("getNCells", &MeshWrapper::getNCells)
  .def("getDimension", &MeshWrapper::getDimension)
	.def(bp::self_ns::str(bp::self_ns::self));
  bp::register_ptr_to_python< std::shared_ptr<MeshWrapper> >();
  bp::implicitly_convertible<std::shared_ptr<MeshWrapper>,std::shared_ptr<mesh::MeshInterface> >();
}

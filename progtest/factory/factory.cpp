#include  <boost/python.hpp>
#include  <boost/python/numpy.hpp>

template<class T> T * get_pointer( std::shared_ptr<T> const& p ) {return p.get();}

class A
{
private:
  int _content;
public:
  A() : _content(11) {}
  A(std::string const& s) : _content(11) {}
  int content() const {return _content;}
  static std::shared_ptr<A> create(std::string const& s);
};

std::shared_ptr<A> A::create(std::string const& s)
{
  return std::make_shared<A>(s);
}

BOOST_PYTHON_MODULE(simfemfactory)
{
    using namespace boost::python;
    class_<A, boost::noncopyable>("A", no_init)
    .def("content", &A::content)
    .def("create", &A::create )
    .staticmethod("create");
    register_ptr_to_python< std::shared_ptr<A> >();
    // implicitly_convertible<std::shared_ptr<A>,std::shared_ptr<const A> >();
}

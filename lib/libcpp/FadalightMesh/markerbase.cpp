#include  "FadalightMesh/markerbase.hpp"
#ifdef  CLANG
#include  <numeric>
#else
#include  <ext/numeric>
#endif
#include  <fstream>
#include  <algorithm>

using namespace FadalightMesh;
using namespace std;

/*--------------------------------------------------------------------------*/

MarkerBase::MarkerBase(double marking_parameter) : MarkerInterface(), _marking_parameter(marking_parameter)
{}

/*--------------------------------------------------------------------------*/

class Compare
{
protected:
  const alat::armavec& _e;
public:
  Compare(const alat::armavec& e) : _e(e) {}
  bool operator()(int i, int j) const
  {
    return _e[j] < _e[i];
  }
};

/*--------------------------------------------------------------------------*/

void MarkerBase::_sort()
{
  int n = _indicator.size();
  _C.set_size(n);
  iota(_C.begin(), _C.end(), 0);
  sort( _C.begin(), _C.end(), Compare(_indicator) );
}

/*--------------------------------------------------------------------------*/

void MarkerBase::write(std::string outfilename, arma::file_type datatype)
{
  std::ofstream file( outfilename.c_str() );
  if( not file.is_open() )
  {
    std::cerr<<"*** ERROR in MarkerBase::write(): cannot open file \""<<outfilename<<"\"\n";
    exit(7);
  }
  _marked_cells.save(file, datatype);
  file.close();
}


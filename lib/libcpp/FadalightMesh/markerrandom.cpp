#include  "FadalightMesh/markerrandom.hpp"
#ifdef  CLANG
#include  <numeric>
#include  <algorithm>
#include  <cstdlib>
#else
#include  <ext/numeric>
#include  <ext/algorithm>
#endif
#include  <fstream>

using namespace FadalightMesh;

/*--------------------------------------------------------------------------*/

void MarkerRandom::_mark()
{
  alat::armaivec all(_n);
#ifdef  CLANG
  std::iota(all.begin(), all.end(), 0);
#else
  std::iota(all.begin(), all.end(), 0);
#endif
  int nmarked = (int) ( _marking_parameter*(double)_n );
  nmarked = std::min(_n,nmarked);
  _marked_cells.set_size(nmarked);
// #ifdef  CLANG
  std::random_shuffle( all.begin(), all.end());
  std::copy(all.begin(),all.begin()+nmarked,_marked_cells.begin());
  // std::cerr << "_marking_parameter = " << _marking_parameter << "\n";
  // std::cerr << "nmarked = " << nmarked << "\n";
  // std::cerr << "all:\n";
  // copy(all.begin(), all.end(), std::ostream_iterator<double>(std::cerr, " "));
  // std::cerr << "\n_marked_cells:\n";
  // copy(_marked_cells.begin(), _marked_cells.end(), std::ostream_iterator<double>(std::cerr, " "));
// #else
//   random_sample( all.begin(), all.end(), _marked_cells.begin(), _marked_cells.end() );
// #endif
}

/*--------------------------------------------------------------------------*/

void MarkerRandom::write(std::string outfilename, arma::file_type datatype)
{
  _mark();
  MarkerBase::write(outfilename, datatype);
}

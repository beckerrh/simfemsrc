#include  "FadalightMesh/measureofcell.hpp"
#include  <fstream>
#include  <cassert>

using namespace FadalightMesh;
using namespace std;

/*--------------------------------------------------------------------------*/
void MeasureOfCell::load(std::string filename)
{
  ifstream file( filename.c_str() );
  if( not file.is_open() )
  {
    cerr << "*** MeasureOfCell::read() : cannot read file \"" << filename << "\"\n";
    assert(0);
  }
  alat::armavec::load(file);
}

/*--------------------------------------------------------------------------*/
void MeasureOfCell::save(std::string filename, arma::file_type datatype) const
{
  ofstream file( filename.c_str() );
  if( not file.is_open() )
  {
    cerr << "*** MeasureOfCell::write() : cannot write file \"" << filename << "\"\n";
    assert(0);
  }
  alat::armavec::save(file, datatype);
}

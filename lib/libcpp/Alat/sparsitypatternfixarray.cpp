#include  "Alat/sparsitypatternfixarray.hpp"
#include  <fstream>

using namespace alat;
using namespace std;

/*-------------------------------------------------------------*/

template<int N>
SparsityPatternFixArray<N>::~SparsityPatternFixArray()
{}

/*-------------------------------------------------------------*/

template<int N>
SparsityPatternFixArray<N>::SparsityPatternFixArray() : _col(), _rowstart() {}

/*-------------------------------------------------------------*/

template<int N>
SparsityPatternFixArray<N>::SparsityPatternFixArray(const SparsityPatternFixArray<N>& S) : _col(S._col), _rowstart(S._rowstart) {}

/*-------------------------------------------------------------*/

template<int N>
SparsityPatternFixArray<N>& SparsityPatternFixArray<N>::operator=(const SparsityPatternFixArray<N>& S)
{
  std::cerr << "COUCOU\n";
  _col.set_size(S._col.size());
  _col = S._col;
  _rowstart.set_size(S._rowstart.size());
  _rowstart = S._rowstart;
  return *this;
}

/*-------------------------------------------------------------*/

template<int N>
void SparsityPatternFixArray<N>::set_size(const SparsityPatternFixArraySoft<N>& SP)
{
  int n = SP.size();
  int ntotal = SP.ntotal();
  _col.set_size(ntotal);
  _rowstart.set_size(n+1);

  _rowstart[0] = 0;
  for(int i = 0; i < n; i++)
  {
    _rowstart[i+1] = _rowstart[i]+SP[i].size();
  }
  int pos = 0;
  for(int i = 0; i < n; i++)
  {
    for(typename SparsityPatternFixArraySoft<N>::const_iterator p = SP[i].begin(); p != SP[i].end(); p++)
    {
      _col[pos++] = *p;
    }
  }
}

/*-------------------------------------------------------------*/
template<int N>
void SparsityPatternFixArray<N>::save( std::ostream& out, arma::file_type datatype) const
{
  col().save(out, datatype);
  rowstart().save(out, datatype);
}

/*-------------------------------------------------------------*/
template<int N>
void SparsityPatternFixArray<N>::load(std::istream& in)
{
  col().load(in);
  rowstart().load(in);
}

// /*-------------------------------------------------------------*/
// template<int N>
// void SparsityPatternFixArray<N>::write(const std::string& filename, std::string datatype) const
// {
//   // std::string filename = basename + ".fadalightconnect";
//   std::ofstream file( filename.c_str() );
//   assert( file.is_open() );
//   save(file, datatype);
//   // col().write(file, datatype);
//   // rowstart().write(file, datatype);
//   file.close();
// }
//
// /*-------------------------------------------------------------*/
// template<int N>
// void SparsityPatternFixArray<N>::read(const std::string& filename)
// {
//   // std::string filename = basename + ".fadalightconnect";
//   std::ifstream file( filename.c_str() );
//   if( !file.is_open() )
//   {
//     std::cerr << "*** ERROR in alat::SparsityPattern::read(): could not open file " << filename << "\n";
//     assert(0);
//     exit(1);
//   }
//   load(file);
//   // col().read(file);
//   // rowstart().read(file);
//   file.close();
// }

/*-------------------------------------------------------------*/
template<int N>
int  SparsityPatternFixArray<N>::n() const
{
  return _rowstart.size()-1;
}

template<int N>
int  SparsityPatternFixArray<N>::ntotal() const
{
  return _col.size();
}

template<int N>
int  SparsityPatternFixArray<N>::rowsize(int i) const
{
  return _rowstart[i+1]-_rowstart[i];
}

template<int N>
const alat::Vector<alat::FixArray<N, int> >&  SparsityPatternFixArray<N>::col() const
{
  return _col;
}

template<int N>
alat::Vector<alat::FixArray<N, int> >&  SparsityPatternFixArray<N>::col()
{
  return _col;
}

template<int N>
const alat::armaivec&  SparsityPatternFixArray<N>::rowstart() const
{
  return _rowstart;
}

template<int N>
alat::armaivec& SparsityPatternFixArray<N>::rowstart()
{
  return _rowstart;
}

template<int N>
alat::FixArray<N, int>&  SparsityPatternFixArray<N>::col(int pos)
{
  return _col[pos];
}

template<int N>
const alat::FixArray<N, int>&  SparsityPatternFixArray<N>::col(int pos) const
{
  return _col[pos];
}

// template<int N>
// int&  SparsityPatternFixArray<N>::rowstart(int i)
// {
//   return _rowstart[i];
// }

template<int N>
int SparsityPatternFixArray<N>::rowstart(int i) const
{
  return _rowstart[i];
}

// template<int N>
// int&  SparsityPatternFixArray<N>::rowstop(int i)
// {
//   return _rowstart[i+1];
// }

template<int N>
int  SparsityPatternFixArray<N>::rowstop(int i) const
{
  return _rowstart[i+1];
}

template<int N>
std::ostream& alat::operator<<(std::ostream& s, const SparsityPatternFixArray<N>& A)
{
  s << "start:\n"<< A.rowstart() << std::endl;
  s << "col:\n"<< A.col() << std::endl;
  return s;
}

/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/

template class alat::SparsityPatternFixArray<2>;

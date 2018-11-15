#include  "Alat/sparsitypattern.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Alat/set.hpp"
#include  <algorithm>
#include  <cassert>
#include  <fstream>

using namespace alat;
using namespace std;

/*-------------------------------------------------------------*/
SparsityPattern::~SparsityPattern() {}
SparsityPattern::SparsityPattern() : _col(), _rowstart() {}
SparsityPattern::SparsityPattern(const SparsityPattern& S) : _col( S.col() ), _rowstart( S.rowstart() ) {}
SparsityPattern::SparsityPattern(const std::vector<aramintvec>& vectors)
{
  int n = vectors.size();
  int ntotal=0;
  for(int i=0;i<n;i++) {ntotal += vectors[i].size();}
  _col.set_size( ntotal );
  _rowstart.set_size(n+1);

  _rowstart[0] = 0;
  for(int i = 0; i < n; i++)
  {
    _rowstart[i+1] = _rowstart[i]+vectors[i].size();
  }
  int pos = 0;
  for(int i = 0; i < n; i++)
  {
    for(int ii = 0; ii < vectors[i].size(); ii++)
    {
      _col[pos++] = vectors[i][ii];
    }
  }
}
SparsityPattern::SparsityPattern(int colsize, int rowstartsize)
{
  _col.set_size(colsize);
  _rowstart.set_size(rowstartsize);
}
SparsityPattern& SparsityPattern::operator=(const SparsityPattern& S)
{
  _rowstart.set_size( S.rowstart().size() );
  _rowstart = S.rowstart();
  _col.set_size( S.col().size() );
  _col = S.col();
  return *this;
}

/*--------------------------------------------------------------------------*/
std::ostream& alat::operator<<(std::ostream& os, const SparsityPattern& A)
{
  os << "start:\n"<< A.rowstart().t() << std::endl;
  os << "col:\n"<< A.col().t() << std::endl;
  for(int i = 0; i < A.n(); i++)
  {
    os << "\nline " << i << " : ";
    for(int pos = A.rowstart(i); pos < A.rowstop(i); pos++)
    {
      os << A.col(pos) << " ";
    }
  }
  return os;
}

int SparsityPattern::n() const
{
  if(_rowstart.size() == 0)
  {
    return 0;
  }
  return _rowstart.size()-1;
}

int SparsityPattern::ntotal() const
{
  return _col.size();
}

int SparsityPattern::rowsize(int i) const
{
  if(_rowstart.size() == 0)
  {
    return 0;
  }
  return _rowstart[i+1]-_rowstart[i];
}

/*-------------------------------------------------------------*/
const alat::SparsityPattern::aramintvec&  SparsityPattern::col() const{return _col;}
alat::SparsityPattern::aramintvec& SparsityPattern::col(){return _col;}
const alat::SparsityPattern::aramintvec&  SparsityPattern::rowstart() const{return _rowstart;}
alat::SparsityPattern::aramintvec& SparsityPattern::rowstart(){return _rowstart;}
int SparsityPattern::col(int pos) const{return _col[pos];}
int SparsityPattern::rowstart(int i) const{return _rowstart[i];}
int  SparsityPattern::rowstop(int i) const{return _rowstart[i+1];}
int SparsityPattern::get(int i, int ii) const{return _col[_rowstart[i]+ii];}

/*-------------------------------------------------------------*/
void SparsityPattern::set_size(int n, int nt)
{
  _col.set_size(nt);
  _rowstart.set_size(n+1);
}

/*-------------------------------------------------------------*/
void SparsityPattern::set_size(const SparsityPattern& SP)
{
  set_size( SP.n(), SP.ntotal() );
  _col = SP.col();
  _rowstart = SP.rowstart();
  sort();
}

/*-------------------------------------------------------------*/
void SparsityPattern::set_size(const SparsityPatternSoft& sparsitypatternsoft)
{
  // std::cerr << "SparsityPattern::set_size()" << sparsitypatternsoft << "\n";
  int n = sparsitypatternsoft.size();
  _col.set_size( sparsitypatternsoft.ntotal() );
  _rowstart.set_size(n+1);
  _rowstart[0] = 0;
  for(int i = 0; i < n; i++){_rowstart[i+1] = _rowstart[i]+sparsitypatternsoft[i].size();}
  int pos = 0;
  for(int i = 0; i < n; i++)
  {
    for(alat::IntSet::const_iterator p = sparsitypatternsoft[i].begin(); p != sparsitypatternsoft[i].end(); p++){_col[pos++] = *p;}
  }
  sort();
}

/*-------------------------------------------------------------*/
void SparsityPattern::set_size(const alat::Vector<aramintvec>& SP)
{
  // std::cerr << "::::::" << SP << "\n";
  int n = SP.size();
  int ntotal = 0;
  for(int i = 0; i < n; i++)
  {
    ntotal += SP[i].size();
  }
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
    for(aramintvec::const_iterator p = SP[i].begin(); p != SP[i].end(); p++)
    {
      _col[pos++] = *p;
    }
  }
}
/*-------------------------------------------------------------*/
void SparsityPattern::print( std::ostream& os) const
{
  os << "SparsityPattern n=" << n() << " ntotal " << ntotal() << "\n";
  for(int i = 0; i < n(); i++)
  {
    for(int pos = rowstart(i); pos < rowstop(i); pos++)
    {
      os << col(pos) << " ";
    }
    os << "\n";
  }
}

/*-------------------------------------------------------------*/
void SparsityPattern::save(std::ostream& out, arma::file_type datatype) const
{
  col().save(out, datatype);
  rowstart().save(out, datatype);
}

/*-------------------------------------------------------------*/
void SparsityPattern::load(std::istream& in)
{
  col().load(in);
  rowstart().load(in);
}

/*-------------------------------------------------------------*/
void SparsityPattern::save(const std::string& filename, arma::file_type datatype) const
{
  // std::string filename = basename + ".fadalightconnect";
  std::ofstream file( filename.c_str() );
  assert( file.is_open() );
  save(file, datatype);
  // col().write(file, datatype);
  // rowstart().write(file, datatype);
  file.close();
}

/*-------------------------------------------------------------*/
void SparsityPattern::load(const std::string& filename)
{
  // std::string filename = basename + ".fadalightconnect";
  std::ifstream file( filename.c_str() );
  if( !file.is_open() )
  {
    std::cerr << "*** ERROR in alat::SparsityPattern::read(): could not open file " << filename << "\n";
    assert(0);
    exit(1);
  }
  load(file);
  // col().read(file);
  // rowstart().read(file);
  file.close();
}

/*--------------------------------------------------------------------------*/
void SparsityPattern::enlarge(int enlarge_stencil)
{
  if(enlarge_stencil == 0)
  {
    return;
  }
  alat::SparsityPatternSoft SPS( n() );
  for(int i = 0; i < n(); i++)
  {
    for(int pos = rowstart(i); pos < rowstop(i); pos++)
    {
      SPS[i].insert( col(pos) );
    }
  }

  for(int ilarge = 0; ilarge < enlarge_stencil; ilarge++)
  {
    alat::SparsityPatternSoft SPSadd( n() );
    for(int i = 0; i < n(); i++)
    {
      for(std::set<int>::const_iterator it = SPS[i].begin(); it != SPS[i].end(); it++)
      {
        int j = *it;
        for(std::set<int>::const_iterator it2 = SPS[j].begin(); it2 != SPS[j].end(); it2++)
        {
          SPSadd[i].insert(*it2);
        }
      }
    }
    for(int i = 0; i < n(); i++)
    {
      for(std::set<int>::const_iterator it = SPSadd[i].begin(); it != SPSadd[i].end(); it++)
      {
        SPS[i].insert(*it);
      }
    }
  }
  set_size(SPS);
}

/*--------------------------------------------------------------------------*/
int SparsityPattern::find(int i, int j) const
{
  for(int pos = rowstart(i); pos < rowstop(i); pos++)
  {
    if(col(pos) == j)
    {
      return pos;
    }
  }
  std::cerr << "*** ERROR alat::SparsityPattern::find()";
  std::cerr << "no such coupling: "<< i <<" "<<j<<std::endl;
  assert(0);
  return -1;
}

/*--------------------------------------------------------------------------*/
void SparsityPattern::sort()
{
  for(int i = 0; i < n(); i++)
  {
    std::sort(&_col[_rowstart[i]], &_col[_rowstart[i+1]]);
  }
}

/*--------------------------------------------------------------------------*/
void SparsityPattern::clear()
{
  _col.clear();
  _rowstart.clear();
}

/*--------------------------------------------------------------------------*/
void SparsityPattern::reconstructWithNumbering(const alat::SparsityPattern* sparsitypattern, const alat::SparsityPattern::aramintvec& p, const alat::SparsityPattern::aramintvec& pinv,  int enlarge_stencil)
{
  assert(sparsitypattern);
  // if(enlarge_stencil)
  int n = sparsitypattern->n();
  alat::SparsityPatternSoft sparsitypatternsoft(n);

  for(int i = 0; i < n; i++)
  {
    int pi = p[i];
    sparsitypatternsoft[i].insert(i);
    for(int pos = sparsitypattern->rowstart(pi); pos < sparsitypattern->rowstop(pi); pos++)
    {
      sparsitypatternsoft[i].insert(pinv[sparsitypattern->col(pos)]);
    }
  }
  // std::cerr << "enlarge_stencil " << enlarge_stencil << "\n";
  for(int ilarge = 0; ilarge < enlarge_stencil; ilarge++)
  {
    alat::SparsityPatternSoft sparsitypatternsoftbig(n);
    for(int i = 0; i < n; i++)
    {
      for(std::set<int>::const_iterator it = sparsitypatternsoft[i].begin(); it != sparsitypatternsoft[i].end(); it++)
      {
        int j = *it;
        for(std::set<int>::const_iterator it2 = sparsitypatternsoft[j].begin(); it2 != sparsitypatternsoft[j].end(); it2++)
        {
          sparsitypatternsoftbig[i].insert(*it2);
        }
      }
    }
    for(int i = 0; i < n; i++)
    {
      for(std::set<int>::const_iterator it = sparsitypatternsoftbig[i].begin(); it != sparsitypatternsoftbig[i].end(); it++)
      {
        sparsitypatternsoft[i].insert(*it);
      }
    }
  }
  set_size(sparsitypatternsoft);
}

/*--------------------------------------------------------------------------*/
void SparsityPattern::setDiagonal(aramintvec& diagonal) const
{
  diagonal.set_size(n());
  for(int i = 0; i < n(); i++)
  {
    bool found = 0;
    for(int pos = rowstart(i); pos < rowstop(i); pos++)
    {
      if(col(pos) == i)
      {
        diagonal[i] = pos;
        found = 1;
        continue;
      }
    }
    if(!found)
    {
      std::cerr << "*** ERROR in IluBlockmatrix::_setDiag() : i = " << i << "\n";
      for(int pos = rowstart(i); pos < rowstop(i); pos++)
      {
        std::cerr  << col(pos) << " ";
      }
      std::cerr << "\n";
      assert(0);
      exit(1);
    }
  }
}

#ifndef  __Alat_SparsityPatternFixArray_h
#define  __Alat_SparsityPatternFixArray_h

#include  "sparsitypatternfixarraysoft.hpp"

/*-----------------------------------------------------------------*/

namespace alat
{
  template<int N>
  class SparsityPatternFixArray
  {
protected:
    alat::armaivec _rowstart;
    alat::Vector<alat::FixArray<N, int> >   _col;

public:
    ~SparsityPatternFixArray();
    SparsityPatternFixArray();
    SparsityPatternFixArray(const SparsityPatternFixArray& S);
    SparsityPatternFixArray& operator=(const SparsityPatternFixArray& S);

    int  n() const;
    int  ntotal() const;
    int  rowsize(int i) const;
    const alat::Vector<alat::FixArray<N, int> >&  col() const;
    alat::Vector<alat::FixArray<N, int> >&  col();
    const alat::armaivec&  rowstart() const;
    alat::armaivec& rowstart();
    alat::FixArray<N, int>&  col(int pos);
    const alat::FixArray<N, int>&  col(int pos) const;
    // int&  rowstart(int i);
    int rowstart(int i) const;
    // int&  rowstop(int i);
    int  rowstop(int i) const;
    void set_size(const SparsityPatternFixArraySoft<N>& SP);

    // void write(const std::string& filename, std::string datatype = "binary") const;
    // void read(const std::string& filename);
    void save( std::ostream& out, arma::file_type datatype = arma::arma_binary) const;
    void load(std::istream& in);
  };
  template<int N>
  std::ostream& operator<<(std::ostream& s, const SparsityPatternFixArray<N>& A);
}

#endif

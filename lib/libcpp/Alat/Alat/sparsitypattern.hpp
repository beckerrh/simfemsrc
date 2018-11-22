#ifndef  __Alat_SparsityPattern_h
#define  __Alat_SparsityPattern_h

#include  "Alat/armadillo.hpp"

/*-----------------------------------------------------------------*/

namespace alat
{
  template<class T>
  class Vector;
  class SparsityPatternSoft;

  class SparsityPattern
  {
  public:
    #ifdef SPARSITYPATTERNLONG
    typedef arma::Col<long> aramintvec;
    #else
    typedef arma::Col<int> aramintvec;
    #endif

  protected:
    aramintvec _col, _rowstart;

  public:
    ~SparsityPattern();
    SparsityPattern();
    SparsityPattern(const SparsityPattern& sparsitypattern);
    SparsityPattern(const std::vector<aramintvec>& vectors);
    SparsityPattern(int colsize, int rowstartsize);
    SparsityPattern& operator=(const SparsityPattern& S);

    void clear();
    int n() const;
    int ntotal() const;
    int rowsize(int i) const;
    const aramintvec&  col() const;
    aramintvec&  col();
    const aramintvec&  rowstart() const;
    aramintvec& rowstart();
    int col(int pos) const;
    int rowstart(int i) const;
    int rowstop(int i) const;
    int get(int i, int ii) const;

    void set_size(int n, int nt);
    void set_size(const SparsityPattern& SP);
    void set_size(const SparsityPatternSoft& SP);
    void set_size(const alat::Vector<aramintvec>& SP);

    int find(int i, int j) const;
    void sort();
    void print( std::ostream& os) const;
    void save(std::ostream& out, arma::file_type datatype = arma::arma_binary) const;
    void load(std::istream& in);
    void save(const std::string& filename, arma::file_type datatype = arma::arma_binary) const;
    void load(const std::string& filename);
    void enlarge(int enlarge_stencil);
    void reconstructWithNumbering(const alat::SparsityPattern* sparsitypattern, const aramintvec& p, const aramintvec& pinv,  int enlarge_stencil = 0);
    void setDiagonal(aramintvec& diagonal) const;
  };
  std::ostream& operator<<(std::ostream& s, const SparsityPattern& A);
}

#endif

#include  "Alat/sparsitypatternfixarraysoft.hpp"

using namespace alat;
using namespace std;

/*-------------------------------------------------------------*/

template<int N>
SparsityPatternFixArraySoft<N>::~SparsityPatternFixArraySoft()
{}

/*-------------------------------------------------------------*/

template<int N>
SparsityPatternFixArraySoft<N>::SparsityPatternFixArraySoft() : alat::Vector<alat::Set<alat::FixArray<N, int> > >()
{}

/*-------------------------------------------------------------*/

template<int N>
SparsityPatternFixArraySoft<N>::SparsityPatternFixArraySoft(const SparsityPatternFixArraySoft& S) : alat::Vector<alat::Set<alat::FixArray<N, int> > >(S)
{}

/*-------------------------------------------------------------*/

template<int N>
SparsityPatternFixArraySoft<N>::SparsityPatternFixArraySoft(int n) : alat::Vector<alat::Set<alat::FixArray<N, int> > >(n) {}

/*-------------------------------------------------------------*/

template<int N>
SparsityPatternFixArraySoft<N>& SparsityPatternFixArraySoft<N>::operator=(const SparsityPatternFixArraySoft<N>& S)
{
  alat::Vector<alat::Set<alat::FixArray<N, int> > >::operator=(S);
  return *this;
}

/*-------------------------------------------------------------*/

template<int N>
std::string SparsityPatternFixArraySoft<N>::getName() const
{
  return "SparsityPatternFixArraySoft";
}

/*-------------------------------------------------------------*/

template<int N>
int SparsityPatternFixArraySoft<N>::ntotal() const
{
  int n = alat::Vector<alat::Set<alat::FixArray<N, int> > >::size();
  int ntotal = 0;
  for(int i = 0; i < n; i++)
  {
    ntotal += ( *this )[i].size();
  }
  return ntotal;
}

/*-------------------------------------------------------------*/

template<int N>
int SparsityPatternFixArraySoft<N>::rowsize(int i) const
{
  return ( *this )[i].size();
}

/*-------------------------------------------------------------*/
/*-------------------------------------------------------------*/

template class alat::SparsityPatternFixArraySoft<2>;

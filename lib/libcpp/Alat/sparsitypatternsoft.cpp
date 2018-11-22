#include  "Alat/sparsitypatternsoft.hpp"
#include  <cassert>

using namespace alat;

/*---------------------------------------------------------*/
SparsityPatternSoft::SparsityPatternSoft() : alat::Vector<alat::IntSet>( ) {}
SparsityPatternSoft::SparsityPatternSoft(int n) : alat::Vector<alat::IntSet>(n) {}
SparsityPatternSoft::SparsityPatternSoft(const SparsityPatternSoft& s) : alat::Vector<alat::IntSet>(s) {}
std::string SparsityPatternSoft::getClassName() const
{
  return "SparsityPatternSoft";
}

/*---------------------------------------------------------*/
int SparsityPatternSoft::ntotal() const
{
  int n = alat::Vector<alat::IntSet>::size();
  int ntotal = 0;
  for(int i = 0; i < n; i++)
  {
    ntotal += rowsize(i);
  }
  return ntotal;
}

/*---------------------------------------------------------*/

int SparsityPatternSoft::rowsize(int i) const
{
  return ( *this )[i].size();
}

void SparsityPatternSoft::set_size(int n)
{
  alat::Vector<alat::IntSet>::set_size(n);
}

/*---------------------------------------------------------*/

void SparsityPatternSoft::reInitTranspose(int n, const SparsityPatternSoft& s)
{
  alat::Vector<alat::IntSet>::set_size(n);
  for(size_t i = 0; i < s.size(); i++)
  {
    for(std::set<int>::const_iterator q = s[i].begin(); q != s[i].end(); q++)
    {
      assert(*q < n);
      alat::Vector<alat::IntSet>::operator[](* q).insert(i);
    }
  }
}

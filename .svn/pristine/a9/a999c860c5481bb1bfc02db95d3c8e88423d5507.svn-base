#ifndef __Alat_SparsityPatternFixArraySoft_h
#define __Alat_SparsityPatternFixArraySoft_h

#include  "Alat/fixarray.hpp"
#include  "Alat/set.hpp"
#include  "Alat/vector.hpp"

/*---------------------------------------------------------*/

namespace alat
{
  template<int N>
  class SparsityPatternFixArraySoft : public alat::Vector<alat::Set<alat::FixArray<N, int> > >
  {
public:
    typedef typename alat::Set<alat::FixArray<N, int> >::const_iterator const_iterator;

public:
    ~SparsityPatternFixArraySoft();
    SparsityPatternFixArraySoft();
    SparsityPatternFixArraySoft(const SparsityPatternFixArraySoft& S);
    SparsityPatternFixArraySoft(int n);
    SparsityPatternFixArraySoft& operator=(const SparsityPatternFixArraySoft& S);

    std::string getName() const;

    int ntotal() const;
    int rowsize(int i) const;
  };
}

/*---------------------------------------------------------*/

#endif

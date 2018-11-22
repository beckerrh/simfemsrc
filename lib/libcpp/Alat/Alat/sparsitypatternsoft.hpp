#ifndef __Alat_SparsityPatternSoft_h
#define __Alat_SparsityPatternSoft_h

#include  "Alat/set.hpp"
#include  "Alat/vector.hpp"

/*---------------------------------------------------------*/

namespace alat
{
  class SparsityPatternSoft : public alat::Vector<alat::IntSet>
  {
public:
    SparsityPatternSoft();
    SparsityPatternSoft(int n);
    SparsityPatternSoft(const SparsityPatternSoft& s);
    std::string getClassName() const;

    int ntotal() const;
    int  rowsize(int i) const;
    void set_size(int n);
    void reInitTranspose(int n, const SparsityPatternSoft& s);
  };
}

/*---------------------------------------------------------*/

#endif

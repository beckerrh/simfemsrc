#ifndef __Alat_FixArray_h
#define __Alat_FixArray_h

#include  <iostream>
#include  <iterator>
#include  <string>
#include  <cassert>
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/

namespace alat
{
  template<int N, class T>
  class FixArray
  {
public:
    typedef  T*        iterator;
    typedef  const T*  const_iterator;

private:
    T val[N];

public:
    ~FixArray<N, T>();
    FixArray<N, T>();
    FixArray<N, T>(const T &d);
    FixArray<N, T>( const FixArray<N, T>& fixarray);
    FixArray<N, T>( const alat::armaivec& vec);
    FixArray<N, T>& operator=( const FixArray<N, T>& fixarray);
    FixArray<N, T>& operator=( const T& t);
    std::string getClassName() const;
    const T* begin() const;
    const T* end() const;
    T* begin();
    T* end();
    size_t   size()            const;
    const T& operator[](int i) const;
    T&  operator[](int i);
    FixArray<N, T>& operator+=(const FixArray<N, T>& v);
    bool operator<(const FixArray<N, T>& v) const;
    bool operator!=(const FixArray<N, T>& v) const;
    bool operator==(const FixArray<N, T>& v) const;

    void sort();
    friend std::ostream& operator<<(std::ostream& s, const FixArray<N, T>& A)
    {
      copy( A.begin(), A.end(), std::ostream_iterator<T>(s, " ") );
      return s;
    }

    friend std::istream& operator>>(std::istream& s, FixArray<N, T>& A)
    {
      typename FixArray<N, T>::iterator p = A.begin();
      while( p != A.end() )
      {
        s >> *p++;
      }
      return s;
    }
    std::ostream& writeBin(std::ostream& out) const;
    std::istream& readBin (std::istream& in);
    std::ostream& save(std::ostream& out, arma::file_type datatype = arma::arma_binary) const;
    std::istream& load(std::istream& in);
  };


/*----------------------------------------------------------*/
  template<int N, class T>
  std::ostream& FixArray<N, T>::writeBin(std::ostream& out) const
  {
    out.write( reinterpret_cast<const char*>( &( *this )[0] ), N*sizeof( T ) );
    return out;
  }
  template<int N, class T>
  std::istream& FixArray<N, T>::readBin(std::istream& in)
  {
    in.read( reinterpret_cast<char*>( &( *this )[0] ), N*sizeof( T ) );
    return in;
  }

/*----------------------------------------------------------*/
  template<int N, class T>
  std::ostream& FixArray<N, T>::save(std::ostream& out, arma::file_type datatype) const
  {
    if(datatype == arma::arma_binary)
    {
      // out<<datatype<<"binary"<<std::endl;
      out<<"binary"<<std::endl;
      writeBin(out);
    }
    else
    {
      // out<<datatype<<"ascii"<<std::endl;
      out<<"ascii"<<std::endl;
      out<<*this;
    }
    return out;
  }
  template<int N, class T>
  std::istream& FixArray<N, T>::load(std::istream& in)
  {
    std::string datatype;
    in>>datatype;
    if(datatype == "ascii")
    {
      in>>*this;
    }
    else if(datatype == "binary")
    {
      std::string str;
      while(str == "")
      {
        getline(in, str);
      }
      readBin(in);
    }
    else
    {
      std::cerr<<"*** FixArray<N, T>::read() : wrong datatype \""<<datatype<<"\"\n";
      assert(0);
    }
    return in;
  }



  /*--------------------------------------------------------------------------*/
  template<int N, class T>
  FixArray<N, T>::~FixArray<N, T>(){}
  template<int N, class T>
  FixArray<N, T>::FixArray(){}
  template<int N, class T>
  FixArray<N, T>::FixArray( const FixArray<N, T>& fixarray )
  {
    std::copy( fixarray.begin(), fixarray.end(), begin() );
  }
  template<int N, class T>
  FixArray<N, T>::FixArray( const alat::armaivec& vec )
  {
    assert(vec.size()==N);
    std::copy( vec.begin(), vec.end(), begin() );
  }
  template<int N, class T>
  FixArray<N, T>::FixArray( const T& d )
  {
    *this = d;
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  FixArray<N, T>& FixArray<N, T>::operator=( const FixArray<N, T>& fixarray)
  {
    std::copy( fixarray.begin(), fixarray.end(), begin() );
    return *this;
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  FixArray<N, T>& FixArray<N, T>::operator=( const T& t)
  {
    iterator p( begin() );
    while( p < end() )
    {
      *p++ = t;
    }
    return *this;
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  std::string FixArray<N, T>::getClassName() const
  {
    return "FixArray";
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  void FixArray<N, T>::sort()
  {
    std::sort(begin(),end());
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  const T* FixArray<N, T>::begin() const
  {
    return &( val[0] );
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  const T* FixArray<N, T>::end() const
  {
    return &( val[0] )+N;
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  T* FixArray<N, T>::begin()
  {
    return &( val[0] );
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  T* FixArray<N, T>::end()
  {
    return &( val[0] )+N;
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  size_t FixArray<N, T>::size() const
  {
    return N;
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  const T& FixArray<N, T>::operator[](int i) const
  {
    return val[i];
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  T& FixArray<N, T>::operator[](int i)
  {
    return val[i];
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  FixArray<N, T>& FixArray<N, T>::operator+=(const FixArray<N, T>& v)
  {
    iterator p( begin() );
    const_iterator q( v.begin() );
    while( p < end() )
    {
      *p++ += *q++;
    }
    return *this;
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  bool FixArray<N, T>::operator<(const FixArray<N, T>& v) const
  {
    const_iterator p( begin() );
    const_iterator q( v.begin() );
    while( p < end() )
    {
      if(*p < *q)
      {
        return true;
      }
      if(*q < *p)
      {
        return false;
      }
      p++;
      q++;
    }
    return false;
  }

  /*--------------------------------------------------------------------------*/

  template<int N, class T>
  bool FixArray<N, T>::operator==(const FixArray<N, T>& v) const
  {
    return not operator!=(v);
  }

  template<int N, class T>
  bool FixArray<N, T>::operator!=(const FixArray<N, T>& v) const
  {
    const_iterator p( begin() );
    const_iterator q( v.begin() );
    while( p < end() )
    {
      if(*p != *q)
      {
        return true;
      }
      p++;
      q++;
    }
    return false;
  }

}

/*--------------------------------------------------------------------------*/

#endif

#ifndef __Alat_List_h
#define __Alat_List_h

#include  <iostream>
#include  <iterator>
#include  <list>

/*----------------------------------------------------------*/

namespace alat
{
  template<class T>
  class List : public std::list<T>
  {
public:
    typedef typename std::list<T>::iterator iterator;
    typedef typename std::list<T>::const_iterator const_iterator;

    ~List();
    List();
    List(const T& v);
    List(const List& v);
    List(const std::list<T>& v);
    template <class InputIterator>
    List(InputIterator first, InputIterator last);
    List& operator=(const List& v);
    friend std::ostream& operator<<(std::ostream& s, const List<T>& A)
    {
      if( !A.size() )
      {
        return s;
      }
      std::ostream_iterator<T>  os(s, " ");
      copy(A.begin(), A.end(), os);
      return s;
    }
  };
  typedef alat::List<int> IntList;

  /*----------------------------------------------------------*/
  template<class T>
  List<T>::~List<T>( )  {}

  /*----------------------------------------------------------*/
  template<class T>
  List<T>::List() : std::list<T>()    {}

  /*----------------------------------------------------------*/
  template<class T>
  List<T>::List(const T& v) : std::list<T>()
  {
    std::list<T>::push_pack(v);
  }

  /*----------------------------------------------------------*/
  template<class T>
  List<T>::List(const List& v) : std::list<T>(v)
  {
    *this = v;
  }

  /*----------------------------------------------------------*/
  template<class T>
  List<T>::List(const std::list<T>& v) : std::list<T>(v)   {}

  /*----------------------------------------------------------*/
  template<class T>
  template <class InputIterator>
  List<T>::List(InputIterator first, InputIterator last) : std::list<T>(first, last)   {}

  /*----------------------------------------------------------*/
  template<class T>
  List<T>& List<T>::operator=(const List<T>& v)
  {
    std::list<T>::operator=(v);
    return *this;
  }
}

#endif

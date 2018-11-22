#include  "Alat/node.hpp"
#include  "Alat/tokenize.hpp"
#include  <cmath>
#include  <iostream>
#include  <stdlib.h>

using namespace alat;

/*--------------------------------------------------------------------------*/
Node::~Node() {}
Node::Node() : arma::vec3()
{
  zeros();
}
Node::Node(const Node& c) : arma::vec3(c) {}
Node::Node(double x, double y, double z) :  arma::vec3()
{
  ( *this )[0] = x;
  ( *this )[1] = y;
  ( *this )[2] = z;
}
Node::Node(std::string valuechain) :  arma::vec3()
{
  reInit(valuechain);
}
Node& Node::operator=(const arma::vec3& c)
{
  arma::vec3::operator=(c);
  return *this;
}
// Node& Node::operator=(double d, const Node& c)
// {
//   arma::vec3::operator=(d,c);
//   return *this;
// }

/*--------------------------------------------------------------------------*/
void Node::reInit(std::string valuechain)
{
  // std::cerr<<"valuechain"<<valuechain<<std::endl;
  alat::StringVector bouts = alat::Tokenize(valuechain, "|");
  if(bouts.size() != 3)
  {
    std::cerr << "*** ERROR in Node::reInit( valuechain = " << valuechain << "\n";
    assert(0);
    exit(1);
  }
  for(int i = 0; i < 3; i++)
  {
    ( *this )[i] = atof( bouts[i].c_str() );
  }
}

/*--------------------------------------------------------------------------*/
double Node::norm() const
{
  return arma::norm(*this);
}

/*--------------------------------------------------------------------------*/
void Node::add(double d, const alat::Node& v)
{
  *this += d*v;
}

/*--------------------------------------------------------------------------*/
void Node::equal(double d1, const alat::Node& v1, double d2, const alat::Node& v2)
{
  x() = d1*v1.x() + d2*v2.x();
  y() = d1*v1.y() + d2*v2.y();
  z() = d1*v1.z() + d2*v2.z();
}

/*--------------------------------------------------------------------------*/
void Node::read(std::istream& is)
{
  is>>( *this );
}
void Node::write(std::ostream& os)
{
  os<<( *this );
}
std::ostream& alat::operator<<(std::ostream& s, const Node& A)
{
  copy( A.begin(), A.end(), std::ostream_iterator<double>(s, " ") );
  return s;
}

std::istream& alat::operator>>(std::istream& s, Node& A)
{
  arma::vec3::iterator p = A.begin();
  while( p != A.end() )
  {
    s >> *p++;
  }
  return s;
}

/*--------------------------------------------------------------------------*/
void Node::scale(double d)
{
  *this *= d;
}

/*--------------------------------------------------------------------------*/
double Node::operator*(const Node& v) const
{
  return arma::dot(*this, v);
}

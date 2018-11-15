#ifndef  __Alat_Node_h
#define  __Alat_Node_h

#include  "Alat/armadillo.hpp"

/*----------------------------------------------------------*/
namespace alat
{
  class Node : public arma::vec3
  {
  public:
    ~Node();
    Node();
    Node(const Node& c);
    Node(double x, double y, double z);
    Node(std::string valuechain);
    Node& operator=(const arma::vec3& c);
    // Node& operator=(double d1, const Node& c);
    // Node& operator=(double d1, const Node& c1, double d2, const Node& c2);
    double operator*(const Node& v) const;
    void scale(double d);

    void add(double d, const alat::Node& v);
    void equal(double d1, const alat::Node& v1, double d2, const alat::Node& v2);
    void reInit(std::string valuechain);
    inline const double& x() const
    {
      return ( *this )[0];
    }
    inline const double& y() const
    {
      return ( *this )[1];
    }
    inline const double& z() const
    {
      return ( *this )[2];
    }
    inline double& x()
    {
      return ( *this )[0];
    }
    inline double& y()
    {
      return ( *this )[1];
    }
    inline double& z()
    {
      return ( *this )[2];
    }
    double norm() const;
    void read(std::istream& is);
    void write(std::ostream& os);
  };
  std::ostream& operator<<(std::ostream& s, const Node& A);
  std::istream& operator>>(std::istream& s, Node& A);
}

#endif

#ifndef __FadalightAdaptiveMesh_Node_h
#define __FadalightAdaptiveMesh_Node_h

#include  "Alat/node.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{
  class Node
  {
  private:

    int _oldid,_id;                
    alat::Node _v;

  public:

    Node();
    Node(const Node& N) : _oldid(N.oldid()), _id(N.id()), _v(N.getNode()) {}
    Node(int id, const alat::Node& v) : _oldid(id), _id(id), _v(v) {}
    Node(std::istream& in);
    std::string getClassName() const {return "Node";}

    const int& id() const {return _id;}
    int& id() {return _id;}
    const int& oldid() const {return _oldid;}
    int& oldid() {return _oldid;}

    const alat::Node& getNode() const {return _v;}
    alat::Node& getNode() {return _v;}
    double x() const {return _v.x();}
    double y() const {return _v.y();}
    double z() const {return _v.z();}

    void write(std::ostream& out, std::string datatype) const;
    void read(std::istream& in, std::string datatype);
    friend std::istream& operator >>(std::istream& in,Node & n); 
  };
}

/*--------------------------------------------------------------------------*/

#endif

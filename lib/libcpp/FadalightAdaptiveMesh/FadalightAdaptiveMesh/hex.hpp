#ifndef __FadalightAdaptiveMesh_Hex_h
#define __FadalightAdaptiveMesh_Hex_h

#include  "typedefs.hpp"
#include  "Alat/fixarray.hpp"
#include  "edge.hpp"
#include  "node.hpp"
#include  "volumeinterface.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{
  class Hex : public VolumeInterface
  {
  private:
    int _oldid, _id, _reftype;
    alat::FixArray<6,face_pointer> _sides;
    alat::FixArray<12,edge_pointer> _edges;
    alat::FixArray<8,Node*> _nodes;

  public:
    ~Hex();
    Hex();
    Hex(int id);
    std::string getClassName() const;
    void set_size(std::istream& in, const std::map<int,Node*>& id2node,const std::map<int,edge_pointer>& id2edge, const std::map<int,face_pointer>& id2face);
    int NumberOfNodes() const;
    int NumberOfFaces() const;
    int NumberOfEdges() const;
    const int& id() const;
    int& id();
    const int& oldid() const;
    int& oldid();
    const int& reftype() const;
    int& reftype();
    const face_pointer& face(int i) const;
    face_pointer& face(int i);
    const edge_pointer& edge(int i) const;
    edge_pointer& edge(int i);
    const Node* node(int i) const;
    Node*& node(int i);
    void write(std::ostream& out, std::string datatype) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

#ifndef __FadalightAdaptiveMesh_Tri_h
#define __FadalightAdaptiveMesh_Tri_h

#include  "Alat/fixarray.hpp"
#include  "Alat/node.hpp"
#include  "typedefs.hpp"
#include  "faceinterface.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{
  class Tri : public FaceInterface
  {
  private:
    int _oldid, _id,_boundaryid;                
    alat::FixArray<3,edge_pointer> _edges;
    alat::FixArray<3,Node*> _nodes;

  public:
    Tri();
    Tri(int id);
    Tri(std::istream& is);
    std::string getClassName() const {return "Tri";}
    void set_size( std::istream& in,const std::map<int,Node*>& id2node, const std::map<int,edge_pointer>& id2edge);
    

    int NumberOfNodes() const {return 3;}
    int NumberOfEdges() const {return 3;}
    int id() const {return _id;}
    int& id() {return _id;}
    const int& oldid() const {return _oldid;}
    int& oldid() {return _oldid;}
    int LocalEdge(const Edge & e);
    
    const edge_pointer& edge(int i) const {return _edges[i];}
    edge_pointer& edge(int i) {return _edges[i];}
    const Node* node(int i) const {return _nodes[i];}
    Node*& node(int i) {return _nodes[i];}
    virtual const int& boundaryid() const {return _boundaryid;}
    virtual int& boundaryid(){return _boundaryid;}

    void write(std::ostream& out, std::string datatype) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

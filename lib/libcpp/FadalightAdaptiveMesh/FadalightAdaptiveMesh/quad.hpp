#ifndef __FadalightAdaptiveMesh_Quad_h
#define __FadalightAdaptiveMesh_Quad_h

#include  "Alat/fixarray.hpp"
#include  "node.hpp"
#include  "typedefs.hpp"
#include  "faceinterface.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{
  class Quad : public FaceInterface
  {
  private:
    int _oldid,_boundaryid;
    alat::FixArray<4,edge_pointer> _edges;
    alat::FixArray<4,Node*> _nodes;
    alat::FixArray<2,volume_pointer> _volumes;

  public:
    Quad();
    Quad(int id);
    Quad(std::istream& is);
    ~Quad(){}
    std::string getClassName() const {return "Quad";}
    void set_size( std::istream& in,const std::map<int,Node*>& id2node, const std::map<int,edge_pointer>& id2edge, const std::map<int, volume_pointer>& id2volume);
    void set_size(std::istream& in, const std::map<int,Node*>& id2node, const std::map<int,edge_pointer>& id2edge);

    int NumberOfNodes() const {return 4;}
    int NumberOfEdges() const {return 4;}
    int id() const {return _id;}
    int& id() {return _id;}
    void setId(int id){_id=id;}
    const int& oldid() const {return _oldid;}
    int& oldid() {return _oldid;}
        
    const edge_pointer& edge(int i) const {return _edges[i];}
    edge_pointer& edge(int i) {return _edges[i];}
    const volume_pointer& volume(int i) const {return _volumes[i];}
    volume_pointer& volume(int i) {return _volumes[i];}
    const Node* node(int i) const {return _nodes[i];}
    Node*& node(int i) {return _nodes[i];}
    const int& boundaryid() const {return _boundaryid;}
    int& boundaryid(){return _boundaryid;}
    void write(std::ostream& out, std::string datatype) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

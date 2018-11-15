#ifndef __FadalightAdaptiveMesh_EdgeRefiner_h
#define __FadalightAdaptiveMesh_EdgeRefiner_h
#include  "FadalightMesh/hanging_info.hpp"
#include "edge.hpp"
/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{ 

//    Nodes
//          
//     1 
//     |
//     .
//     .
//     . 
//     |  
//     3     
//     |
//     2  
//     |   
//     0 
//          
  class EdgeRefiner 
  {
  private: 
    std::vector<Node*> newnodes;
    std::vector<Edge*> newedges;

  public:

    EdgeRefiner(){}
    void  set_size(int n)
    { 
      newnodes.reserve(n-1); newnodes.resize(n-1);
      newedges.reserve(n);newedges.resize(n);
    }
    std::vector<Node*>& getNewNodes(){return newnodes;}
    std::vector<Edge*>& getNewEdges(){return newedges;}     
    void  refine(edge_pointer ep, int& nodecount, int & edgecount);
  };
}

/*--------------------------------------------------------------------------*/

#endif

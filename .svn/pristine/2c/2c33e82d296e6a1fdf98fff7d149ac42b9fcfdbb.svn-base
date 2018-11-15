#ifndef  __FadalightAdaptiveMesh_FaceRefinerInterface_h
#define  __FadalightAdaptiveMesh_FaceRefinerInterface_h

#include  "typedefs.hpp"
#include  "Alat/node.hpp"
#include  <algorithm>

/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{   
 class FaceRefinerInterface 
 {
   protected:
     std::vector<Node*> newnodes;
     std::vector<edge_pointer> newedges;
     std::vector<FaceInterface*> newfaces;
   
   protected:
     FaceRefinerInterface() {}
     
   public:
     std::vector<Node*>& getNewNodes(){return newnodes;}
     std::vector<edge_pointer>& getNewEdges(){return newedges;}
     std::vector<FaceInterface*>& getNewFaces(){return newfaces;}  
     virtual bool  refine(FaceInterface* f,tree<Edge*>& _Edges, int& edgecount){assert(0);}
     virtual bool  coarse(){assert(0);}
     Node correctInternalNodePosition (alat::Node& ninit);
 };
}

/*--------------------------------------------------------------------------*/

#endif

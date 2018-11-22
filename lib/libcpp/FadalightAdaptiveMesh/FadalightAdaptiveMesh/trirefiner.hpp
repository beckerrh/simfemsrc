#ifndef  __FadalightAdaptiveMesh_TriRefiner_h
#define  __FadalightAdaptiveMesh_TriRefiner_h

#include  "typedefs.hpp"
#include  "tri.hpp"
#include  "facerefinerinterface.hpp"
#include  <algorithm>

/*--------------------------------------------------------------------------*/
//
//
//   o
//   | \  
//   |  \  
//   o---o
//   |\  /\  
//   | \/  \ 
//   o--o---o
//
//
//
// 

namespace FadalightAdaptiveMesh
{ 
  class TriRefiner : public FaceRefinerInterface
  {
  public:

    TriRefiner() : FaceRefinerInterface() {}
    void reInit()
    { 
      newedges.reserve(3); newedges.resize(3);
      newfaces.reserve(4); newfaces.resize(4);
    }
    bool refine(FaceInterface* f,tree<Edge*>& _Edges, int& edgecount);
  };
}

/*--------------------------------------------------------------------------*/

#endif

#ifndef __FadalightAdaptiveMesh_TypeDefs_h
#define __FadalightAdaptiveMesh_TypeDefs_h

#include  <map>
#include  <string>
#include  "tree/tree.hh"

/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{
  class VolumeInterface;
  class FaceInterface;
  class Edge;

  typedef tree<VolumeInterface*>::iterator volume_pointer;
  typedef tree<FaceInterface*>::iterator face_pointer;
  typedef tree<Edge*>::iterator edge_pointer;
  typedef tree<VolumeInterface*>::leaf_iterator volume_leafpointer;
  typedef tree<FaceInterface*>::leaf_iterator face_leafpointer;
  typedef tree<Edge*>::leaf_iterator edge_leafpointer;
}

/*--------------------------------------------------------------------------*/

#endif

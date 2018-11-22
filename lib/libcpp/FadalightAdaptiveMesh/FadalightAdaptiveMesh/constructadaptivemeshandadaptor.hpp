#ifndef __FadalightAdaptiveMesh_ConstructAdaptiveMeshAndAdaptor_h
#define __FadalightAdaptiveMesh_ConstructAdaptiveMeshAndAdaptor_h

#include  <string>
#include  "adaptivemeshinterface.hpp"
#include  "meshadaptorinterface.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{
  class ConstructAdaptiveMeshAndAdaptor
  {
  public:
    ConstructAdaptiveMeshAndAdaptor(AdaptiveMeshInterface*& M, MeshAdaptorInterface*& AM, const std::string& meshname, const std::string& adaption_type="refine");
    std::string getClassName() const {return "ConstructAdaptiveMeshAndAdaptor";}
  };
}

/*--------------------------------------------------------------------------*/

#endif

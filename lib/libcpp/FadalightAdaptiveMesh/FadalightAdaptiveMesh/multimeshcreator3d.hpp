#ifndef __FadalightAdaptiveMesh_MultiMeshCreator3d_h
#define __FadalightAdaptiveMesh_MultiMeshCreator3d_h

#include  "FadalightAdaptiveMesh/coarsener3d.hpp"

/*--------------------------------------------------------------------------*/
namespace FadalightAdaptiveMesh
{
  class MultiMeshCreator3d
  {
protected:

   void _coarsen();
   Coarsener3d* _coarsener;
   /// fonction d'accés au pointer sur le maillage adaptatif
   AdaptiveMeshInterface* _getMesh(){return _coarsener->getMesh();}
   /// fonction d'accés à l'arbre des faces
   tree<VolumeInterface*> & _getVolumes(){return _coarsener->getVolumes();}
   void _constructRefineInfo();
public:
    MultiMeshCreator3d(Coarsener3d* c): _coarsener(c){}

    std::string getClassName() const
    {
      return "MultiMeshCreator3d";
    }
   void createMultiMesh(std::string dirname, int nlevels, int ncells, arma::file_type datatype = arma::arma_binary);
   

  };
}

/*--------------------------------------------------------------------------*/

#endif

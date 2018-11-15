#ifndef __FadalightAdaptiveMesh_MultiMeshCreator_h
#define __FadalightAdaptiveMesh_MultiMeshCreator_h

#include "coarsener2d.hpp"

/*--------------------------------------------------------------------------*/
namespace FadalightAdaptiveMesh
{
  class MultiMeshCreator
  {
protected:
   void _coarsen();
   Coarsener2d* _coarsener;
   /// fonction d'accés au pointer sur le maillage adaptatif
   AdaptiveMeshInterface* _getMesh();
   /// fonction d'accés à l'arbre des faces
   tree<FaceInterface*> & _getFaces();
   void _constructRefineInfo();

public:
  ~MultiMeshCreator();
  MultiMeshCreator(Coarsener2d* coarsener);
  MultiMeshCreator(const MultiMeshCreator& multimeshcreator);
  MultiMeshCreator& operator= (const MultiMeshCreator & multimeshcreator);
  std::string getClassName() const;
  void createMultiMesh(std::string dirname, int nlevels, int ncells, arma::file_type datatype = arma::arma_binary);   
  };
}

/*--------------------------------------------------------------------------*/

#endif

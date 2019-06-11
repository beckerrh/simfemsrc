#ifndef __FadalightAdaptiveMesh_BasicAdaptor_h
#define __FadalightAdaptiveMesh_BasicAdaptor_h

#include  "adaptivemeshinterface.hpp"
#include  "meshadaptorinterface.hpp"
#include  "edgerefiner.hpp"
#include  "facerefinerinterface.hpp"
#include  "hexrefiner.hpp"
/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{
  class BasicAdaptor: public MeshAdaptorInterface
  {
protected:
/// pointeur sur le maillage à adapter
   AdaptiveMeshInterface* _adaptive_mesh;
/// outil de raffinement des edges
   EdgeRefiner _edgerefiner;
/// outil de raffinement des faces de type abstrait à spécialiser en fonction des types de faces
   FaceRefinerInterface* _facerefiner;
/// outil de raffinement des faces de type abstrait à spécialiser en fonction des types de faces
   HexRefiner* _volumerefiner;

   /// fonction d'accés à l'ensemble des nodes
   NodeSet & _getNodes(){return _adaptive_mesh->getNodes();}
   /// fonction d'accés à l'arbre des edges
   tree<Edge*> & _getEdges(){return _adaptive_mesh->getEdges();}
   /// fonction d'accés à l'arbre des faces
   tree<FaceInterface*> & _getFaces(){return _adaptive_mesh->getFaces();}
   tree<VolumeInterface*> & _getVolumes(){return _adaptive_mesh->getVolumes();}
   /// fonction d'accés aux informations sur les bords courbes
   // const FadalightMesh::CurvedBoundaryInformation* _getCurvedBoundaries(){return _adaptive_mesh->getCurvedBoundaries();}
   /// fonction d'accès à la dernière identité de node affectée
   int& _getLastNodeId(){return _adaptive_mesh->getLastNodeId();}
   /// fonction d'accès à la dernière identité de node affectée
   int& _getLastEdgeId(){return _adaptive_mesh->getLastEdgeId();}
   /// fonction d'accès à la dernière identité de face affectée
   int& _getLastFaceId(){return _adaptive_mesh->getLastFaceId();}
   /// fonction d'accès à la dernière identité de face affectée
   int& _getLastVolumeId(){return _adaptive_mesh->getLastVolumeId();}
   /// fonction d'accès à la numérotation des edges
   alat::IntMap& _getEdgeId2Id(){return _adaptive_mesh->getEdgeId2Id();}
   /// fonction d'accès à la numérotation des faces
   alat::IntMap& _getFaceId2Id() {return _adaptive_mesh->getFaceId2Id();}
   /// fonction d'accès au cell_map du maillage adaptatif en 2D
   std::map<int, face_pointer>&  _getCellMap2d() {return _adaptive_mesh->getCellMap2d();}
   /// fonction d'accès au cell_map du maillage adaptatif  en 3D
   std::map<int, volume_pointer>&  _getCellMap3d() {return _adaptive_mesh->getCellMap3d();}
   /// fonction d'accès à l'indicateur de mise à jour du cell_map
   bool&  _getCellMapOk(){return _adaptive_mesh->getCellMapOk();}
   /// fonction d'accès à l'indicateur de mise à jour des numérotations
   bool&  _getNumberingOk(){return _adaptive_mesh->getNumberingOk();}

public:
   /// constructeur : on spécialise ici l'outil de refinement des faces en fonction
   /// du type de maillage
    BasicAdaptor(AdaptiveMeshInterface* mesh): _adaptive_mesh(mesh){}
    std::string getClassName() const
    {
      return "BasicAdaptor";
    }
    /// fonction d'accés au pointer sur le maillage adaptatif
    AdaptiveMeshInterface* _getMesh(){return _adaptive_mesh;}
    /// fonction d'accès à la numérotation des nodes
    alat::IntMap& _getNodeId2Id() {return _adaptive_mesh->getNodeId2Id();}
    const alat::IntMap& getEdgeId2Id() const {return _adaptive_mesh->getEdgeId2Id();}
    const alat::IntMap& getFaceId2Id() const {return _adaptive_mesh->getFaceId2Id();}
    const alat::IntMap& getVolumeId2Id() const {return _adaptive_mesh->getVolumeId2Id();}
    tree<FaceInterface*> & getFaces(){return _adaptive_mesh->getFaces();}
    tree<VolumeInterface*> & getVolumes(){return _adaptive_mesh->getVolumes();}
  };
}

/*--------------------------------------------------------------------------*/

#endif

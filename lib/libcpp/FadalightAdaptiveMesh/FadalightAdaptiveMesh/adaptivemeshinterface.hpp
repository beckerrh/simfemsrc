#ifndef __FadalightAdaptiveMesh_AdaptiveMeshInterface_h
#define __FadalightAdaptiveMesh_AdaptiveMeshInterface_h

#include  "FadalightMesh/geometryobject.hpp"
#include  "FadalightMesh/curvedboundaryinformation.hpp"
#include  "faceinterface.hpp"
#include  "Alat/node.hpp"
#include  "edge.hpp"
#include  "typedefs.hpp"
#include  "setdefs.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightAdaptiveMesh
{
  class AdaptiveMeshInterface: public alat::InterfaceBase
  {

public:

    std::string getInterfaceName() const
    {
      return "AdaptiveMeshInterface";
    }
    virtual void readAdaptiveMesh(std::string name, std::string last_in)=0;
    virtual void writeAdaptiveMesh(std::string name, arma::file_type datatype = arma::arma_binary)=0;
    virtual void readFadalightMeshAndInitTrees(std::string name)=0;
    virtual void writeFadalightMesh(std::string name, arma::file_type datatype = arma::arma_binary)=0;
    virtual void constructNumbering()=0;
    virtual void constructCellMap()=0;
    virtual void updateHangingInfo()=0;
    virtual void reInitFadalightMesh()=0;
    virtual NodeSet & getNodes()=0;
    virtual tree<Edge*> & getEdges()=0;
    virtual tree<FaceInterface*> & getFaces()=0;
    virtual tree<VolumeInterface*> & getVolumes()=0;
    virtual const FadalightMesh::CurvedBoundaryInformation* getCurvedBoundaries()=0;
    virtual int& getLastNodeId()=0;
    virtual int& getLastEdgeId()=0;
    virtual int& getLastFaceId()=0;
    virtual int& getLastVolumeId()=0;
    virtual alat::IntMap& getNodeId2Id()=0;
    virtual alat::IntMap& getEdgeId2Id()=0;
    virtual alat::IntMap& getFaceId2Id()=0;
    virtual alat::IntMap& getVolumeId2Id()=0;
    virtual std::map<int, face_pointer>&  getCellMap2d() =0;
    virtual std::map<int, volume_pointer>&  getCellMap3d() =0;
    virtual bool&  getCellMapOk()=0;
    virtual bool&  getNumberingOk()=0;
    virtual FadalightMesh::GeometryObject* getGeometryObject(std::string name) =0;
    virtual bool geometryObjectExists(std::string name) const =0;
    virtual void createGeometryObject(std::string name) =0;
    virtual void writeVtk(std::string name) =0;
  };
}

/*--------------------------------------------------------------------------*/

#endif

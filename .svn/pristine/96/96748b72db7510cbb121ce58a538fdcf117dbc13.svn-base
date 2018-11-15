#ifndef __FadalightMesh_MeshCompositionInterface_h
#define __FadalightMesh_MeshCompositionInterface_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class CouplingMeshInterface;
  class MeshInterface;
  class ModelManagerInterface;

  class MeshCompositionInterface : public alat::InterfaceBase
  {
private:
protected:
    std::string getInterfaceName() const;

public:
    ~MeshCompositionInterface();
    MeshCompositionInterface();
    MeshCompositionInterface( const MeshCompositionInterface& meshcompositioninterface);
    MeshCompositionInterface& operator=( const MeshCompositionInterface& meshcompositioninterface);
    std::string getClassName() const;
    MeshCompositionInterface* clone() const;

    virtual int getDimension() const=0;
    virtual std::string getInfo() const = 0;
    virtual void writeMeshInfo(std::string filename, std::string blockfilename) const;
    // virtual void initCouplingGrids(const alat::armavecModelManagerInterface* _modelmanager);
    virtual void constructFadalightMesh(const std::string& meshname);
    virtual void read(const std::string& basefilename) = 0;
    virtual void write(const std::string& basefilename, arma::file_type datatype = arma::arma_binary) const = 0;
    virtual int getNDomains() const=0;
    virtual int getNCouplingMeshes() const=0;
    virtual const FadalightMesh::MeshInterface* getMesh(int i) const=0;
    virtual void writeH5(std::string filename) const=0;
    virtual const FadalightMesh::CouplingMeshInterface* getCouplingMesh(int i) const=0;
    virtual const FadalightMesh::MeshInterface* getMacroMesh() const;
    virtual int getNCells() const=0;

  };
}

/*--------------------------------------------------------------------------*/

#endif

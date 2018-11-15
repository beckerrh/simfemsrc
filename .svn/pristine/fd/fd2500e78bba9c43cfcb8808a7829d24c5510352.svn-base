#ifndef __FadalightMesh_SingleMeshComposition_h
#define __FadalightMesh_SingleMeshComposition_h

#include  "meshcompositioninterface.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class MeshInterface;
  class SingleMeshComposition : public FadalightMesh::MeshCompositionInterface
  {
private:
    FadalightMesh::MeshInterface* _mesh;

protected:
public:
    ~SingleMeshComposition();
    SingleMeshComposition();
    SingleMeshComposition( const SingleMeshComposition& singlemeshcomposition);
    SingleMeshComposition& operator=( const SingleMeshComposition& singlemeshcomposition);
    std::string getClassName() const;
    SingleMeshComposition* clone() const;

    int getDimension() const;
    std::string getInfo() const;
    void read(const std::string& basefilename);
    void write(const std::string& basefilename, arma::file_type datatype = arma::arma_binary) const;
    void constructFadalightMesh(const std::string& meshname);
    int getNDomains() const;
    int getNCouplingMeshes() const;
    const FadalightMesh::MeshInterface* getMesh(int i) const;
    const FadalightMesh::CouplingMeshInterface* getCouplingMesh(int i) const;
    void writeMeshInfo(std::string filename, std::string blockfilename) const;
    void writeH5(std::string filename) const;
    int getNCells() const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

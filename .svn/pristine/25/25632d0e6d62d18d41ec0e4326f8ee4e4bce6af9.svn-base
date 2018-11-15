#ifndef __FadalightMesh_CouplingMeshInterface_h
#define __FadalightMesh_CouplingMeshInterface_h

#include  "FadalightMesh/meshinterface.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class CouplingMeshInterface : public MeshInterface
  {
private:
protected:
public:
    ~CouplingMeshInterface();
    CouplingMeshInterface();
    CouplingMeshInterface( const CouplingMeshInterface& couplingmeshinterface);
    CouplingMeshInterface& operator=( const CouplingMeshInterface& couplingmeshinterface);
    std::string getClassName() const;
    CouplingMeshInterface* clone() const;

    // virtual void setNLevels(int nlevels) const = 0;
    virtual void computeCouplingSideInformation(int nlevels) const = 0;
  };
}

/*--------------------------------------------------------------------------*/

#endif

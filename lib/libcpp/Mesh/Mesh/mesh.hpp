
#ifndef __Mesh_Mesh_h
#define __Mesh_Mesh_h

#include  "Mesh/meshinterface.hpp"
#include  "Mesh/boundarymeshunit.hpp"
#include  "Mesh/interfacemeshunit.hpp"
#include  "Mesh/interfacemeshunitsmap.hpp"
#include  "Alat/armadillo.hpp"

/*---------------------------------------------------------*/
namespace mesh
{
  class Mesh : public virtual mesh::MeshInterface
  {
  private:
    Mesh();

  protected:
    mesh::MeshUnit _plainmesh;
    int _partion_id;
    bool _construct_bdrymeshes;
    BoundaryMeshUnitsMap _boundarymeshunits;
    InterfaceMeshUnitsMap _interfacemeshunits;

    void constructBdrymesh(BoundaryMeshUnit& mesh, const BoundaryInformation& boundaryinformation, const alat::Map<int, alat::armaimat>& color_to_bdry2);
    void constructInterfacemesh(InterfaceMeshUnit& mesh, const BoundaryInformation& boundaryinformation, const alat::Map<int, alat::armaimat>& color_to_bdry1);
    void _callSidesConstructor(MeshUnitInterface* mesh, const MeshUnitInterface::BoundarySideToColorForSidesConstructor& bsides, const alat::Map<int, alat::armaivec>* partition_to_cells=NULL);

  public:
    ~Mesh();
    Mesh(const Mesh& mesh);
    static std::unique_ptr<MeshInterface> create(meshEnums::meshtype type, int partion_id=1, bool construct_bdrymeshes=false);
    static std::unique_ptr<MeshInterface> create(std::string name, int partion_id=1, bool construct_bdrymeshes=false);
    Mesh& operator=(const Mesh& meshbase);
    void init(std::shared_ptr<const mesh::MeshVisitorInterface> visitor, int partion_id, bool construct_bdrymeshes=false);

    std::string getClassName() const;
    int getPartionId() const;
    void setPartionId(int id);
    void addGeometryObject(meshEnums::geomobjtype type, const GeometryConstructorInterface* geometryconstructor=&trivialgeometryconstructor);
    void addGeometryObjectByName(std::string name);
    std::string getInfo() const;
    int getNCells() const;
    int getDimension() const;

    const MeshUnitInterface* getPlainMesh() const;
    const MeshUnitInterface* getBoundaryMesh(int color) const;
    const MeshUnitInterface* getInterfaceMesh(int color) const;
    const BoundaryMeshUnitsMap& getBoundaryMeshUnitsMap() const;
    const InterfaceMeshUnitsMap& getInterfaceMeshUnitsMap() const;
    MeshUnitInterface* getPlainMesh();
    BoundaryMeshUnitsMap& getBoundaryMeshUnitsMap();
    InterfaceMeshUnitsMap& getInterfaceMeshUnitsMap();

    void saveH5(std::string filename) const;
    void loadH5(std::string filename);
    void writeVtk(std::string filename) const;
    void writeBoundaryVtk(std::string filename) const;
    void readGmsh(std::string filename);
    void exchangeInterfaces(int nproc);
  };
  std::ostream& operator<<(std::ostream &os, const Mesh& mesh);
}

/*---------------------------------------------------------*/

#endif

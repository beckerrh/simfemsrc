#ifndef __FadalightMesh_CoarsenInfo_h
#define __FadalightMesh_CoarsenInfo_h

#include  "FadalightMesh/geometryobject.hpp"
#include  "Alat/sparsitypattern.hpp"

/*---------------------------------------------------------*/

namespace FadalightMesh
{
  /// c'est que pour les quads en ce moments !
  /// cellid -->
  class CoarsenInfo : public FadalightMesh::GeometryObject
  {
protected:
    alat::SparsityPattern _oldnodeids, _oldcellids, _oldsideids, _oldedgeids;
    alat::armaivec _nodenewtoold;

public:
  ~CoarsenInfo();
    CoarsenInfo();
    CoarsenInfo(const CoarsenInfo& coarseninfo);
    CoarsenInfo& operator=( const CoarsenInfo& coarseninfo);
    std::string getClassName() const;

    alat::SparsityPattern& getOldNodes();
    alat::SparsityPattern& getOldCells();
    alat::SparsityPattern& getOldSides();
    alat::SparsityPattern& getOldEdges();
    int getNNodes() const;
    int getNOldNodes(int i) const;
    int getOldNode(int i, int ii) const;
    int getNCells() const;
    int getNOldCells(int i) const;
    int getOldCell(int i, int ii) const;
    int getNSides() const;
    int getNOldSides(int i) const;
    int getOldSide(int i, int ii) const;
    int getNEdges() const;
    int getNOldEdges(int i) const;
    int getOldEdge(int i, int ii) const;
    alat::armaivec& getNodeNewToOld();
    const alat::armaivec& getNodeNewToOld() const;
    void load(std::string filename);
    void save(std::string filename, arma::file_type datatype = arma::arma_binary) const;
  };
}
#endif

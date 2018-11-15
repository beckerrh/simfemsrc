#ifndef __FadalightMesh_RefineInfo_h
#define __FadalightMesh_RefineInfo_h

#include  "FadalightMesh/geometryobject.hpp"
#include  "Alat/sparsitypattern.hpp"
#include  "Alat/fixarray.hpp"

/*---------------------------------------------------------*/

namespace FadalightMesh
{
  /// ce n'est que pour les quads en ce moment !
  class RefineInfo : public FadalightMesh::GeometryObject
  {
protected:
    alat::SparsityPattern _coarsenodeids, _coarsecellids, _coarsesideids, _coarseedgeids;
    alat::armaivec _nodeids, _sideids, _cellids;

public:
    ~RefineInfo();
    RefineInfo();
    RefineInfo(const RefineInfo& refineinfo);
    RefineInfo& operator=( const RefineInfo& refineinfo);
    std::string getClassName() const;

    alat::FixArray<4,int>  refinfoinfonode, refinfoinfoside, refinfoinfocell;

    alat::SparsityPattern& getCoarseNodes();
    alat::SparsityPattern& getCoarseCells();
    alat::SparsityPattern& getCoarseSides();
    alat::SparsityPattern& getCoarseEdges();
    alat::armaivec& getNodeIds();
    alat::armaivec& getSideIds();
    alat::armaivec& getCellIds();
    const alat::SparsityPattern& getCoarseNodes() const;
    const alat::SparsityPattern& getCoarseCells() const;
    const alat::SparsityPattern& getCoarseSides() const;
    const alat::SparsityPattern& getCoarseEdges() const;
    const alat::armaivec& getNodeIds() const;
    const alat::armaivec& getSideIds() const;
    const alat::armaivec& getCellIds() const;
    int getNCoarseCells(int i) const;
    int getNCoarseSides(int i) const;
    int getNCoarseNodes(int i) const;
    int getNCoarseEdges(int i) const;
    int getCoarseCellNumber(int i, int ii) const;
    int getCoarseSideNumber(int i, int ii) const;
    int getCoarseNodesNumber(int i, int ii) const;
    int getCoarseEdgesNumber(int i, int ii) const;
    void load(std::string filename);
    void save(std::string filename, arma::file_type datatype = arma::arma_binary) const;
  };
}
#endif

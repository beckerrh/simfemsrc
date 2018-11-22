#ifndef  __FadalightMesh_CurvedBoundaryInformation_h
#define  __FadalightMesh_CurvedBoundaryInformation_h

#include  "FadalightMesh/curvedboundarydescriptioninterface.hpp"
// #include  "FadalightMesh/curvedboundaryinformationinterface.hpp"
#include  "FadalightMesh/geometryobject.hpp"
#include  "FadalightMesh/meshinterface.hpp"
#include  "Alat/sparsitypatternfixarray.hpp"
#include  "Alat/map.hpp"
#include  <fstream>

/*---------------------------------------------------*/

namespace FadalightMesh
{
  class CurvedBoundaryInformation : public FadalightMesh::GeometryObject
  // , public FadalightMesh::CurvedBoundaryInformationInterface
  {
public:
    typedef alat::Map<int, FadalightMesh::CurvedBoundaryDescriptionInterface*>::const_iterator const_iterator;
    typedef alat::Map<int, FadalightMesh::CurvedBoundaryDescriptionInterface*>::iterator iterator;

private:
    alat::Map<int, FadalightMesh::CurvedBoundaryDescriptionInterface*> _curvedboundaries;
    alat::armaivec _celliscurved;
    alat::SparsityPatternFixArray<2> _spfa;

protected:
  alat::SparsityPatternFixArray<2>& curvedCellIndices();

public:
    ~CurvedBoundaryInformation();
    CurvedBoundaryInformation();
    CurvedBoundaryInformation( const CurvedBoundaryInformation& curvedboundaryinformation);
    CurvedBoundaryInformation& operator=( const CurvedBoundaryInformation& curvedboundaryinformation);
    std::string getClassName() const;

    const alat::Map<int, FadalightMesh::CurvedBoundaryDescriptionInterface*>& get() const;
    const FadalightMesh::CurvedBoundaryDescriptionInterface* get(int color) const;
    FadalightMesh::CurvedBoundaryDescriptionInterface* get(int color);
    FadalightMesh::CurvedBoundaryDescriptionInterface*& getPointer(int color);
    bool boundaryColorIsCurved(int color) const;
    void readCurvedBoundaryDescription(std::istream& in);
    void writeCurvedBoundaryDescription(std::ostream& out) const;
    void load(std::string filename);
    void save(std::string filename, arma::file_type datatype = arma::arma_binary) const;
    void set_size(const CurvedBoundaryInformation& curvedboundaryinformation);

    const alat::SparsityPatternFixArray<2>& curvedCellIndices() const;
    bool cellIsCurved(int iK) const;
    int rowstart(int iK) const;
    int rowstop(int iK) const;
    const alat::FixArray<2, int>& getCurvedInfo(int pos) const;
    void constructBoundaryInformation(const FadalightMesh::MeshInterface* M);
  };
}

/*---------------------------------------------------*/

#endif

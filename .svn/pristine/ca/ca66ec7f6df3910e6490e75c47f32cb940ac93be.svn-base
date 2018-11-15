#ifndef __FadalightMesh_QuadToTri_h
#define __FadalightMesh_QuadToTri_h

#include  "quadrilateralmesh.hpp"
#include  "trianglemesh.hpp"

/*--------------------------------------------------------------------------*/

namespace FadalightMesh
{
  class QuadrilateralMesh;
  class PatchInfo;

  class TriOfQuad
  {
protected:
    alat::Vector<alat::Vector<alat::FixArray<2, int> > > _data;

public:
    void write(std::ostream& out, arma::file_type datatype = arma::arma_binary) const;
    void read(std::istream& in);
    int getNTrianglesOfQuad(int iquad) const;
    int getTriIdOfQuad(int iquad, int iit) const;
    int getLocalIdOfQuadSide(int iquad, int iit) const;
    alat::Vector<alat::Vector<alat::FixArray<2, int> > >& getData();
  };

  /*--------------------------------------------------------------------------*/

  class QuadToTri : public TriangleMesh
  {
protected:
    TriOfQuad _tri_of_quad;
    alat::armaivec _quad_of_tri;
    std::string _infilenamequad;
    QuadrilateralMesh _quadmesh;
    void _makeCrissCross(const FadalightMesh::MeshInterface& QM);

public:
    ~QuadToTri();
    QuadToTri();
    QuadToTri(std::string quadfilename);
    QuadToTri& operator=(const QuadToTri& quadtotri);
    std::string getClassName() const;
    std::string getFileNameQuad() const;
    void convertMesh(const FadalightMesh::MeshInterface& QM, std::string type="crisscross");
    void writeFadalightMesh(const std::string& basefilename, arma::file_type datatype = arma::arma_binary) const;
    void readFadalightMesh(const std::string& basefilename);
    const TriOfQuad& getTriOfQuad() const;
    TriOfQuad& getTriOfQuad();
    alat::Vector<alat::Vector<alat::FixArray<2, int> > >& getTriOfQuadData();
    const alat::armaivec& getQuadOfTri() const;
    const QuadrilateralMesh& getQuadrilateralMesh() const;
    int getCenterNodeIdOfQuad(int iK) const;
    void constructPatchInfo(PatchInfo& patchinfo) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

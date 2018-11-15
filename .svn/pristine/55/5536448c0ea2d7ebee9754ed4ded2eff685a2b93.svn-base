#ifndef __FadalightMesh_HangingSideInfo_h
#define __FadalightMesh_HangingSideInfo_h

#include  "FadalightMesh/geometryobject.hpp"
#include  "Alat/sparsitypattern.hpp"

/*---------------------------------------------------------*/
namespace FadalightMesh
{
  class HangingSideInfo : public FadalightMesh::GeometryObject, public alat::SparsityPattern
  {
  public:
    ~HangingSideInfo();
    HangingSideInfo();
    HangingSideInfo(const HangingSideInfo& hangingsideinfo);
    HangingSideInfo& operator=(const HangingSideInfo& hangingsideinfo);

    std::string getClassName() const;
    int getCellNumber(int i) const;
    int getLocalSide(int i) const;
    int getNumberOfHangingSides(int i) const;
    int getHangingSides(int i, int in) const;
    int& getCellNumber(int i);
    int& getLocalSide(int i);
    int& getHangingSides(int i, int in);
    void set_size(int n_hanging, int n_local_data);
    std::vector<int> getNeighbourCellsOfSide(int iK, int iis, const FadalightMesh::MeshInterface* mesh);
    void load(std::string filename );
    void save(std::string filename, arma::file_type datatype = arma::arma_binary) const;
  };
}


#endif

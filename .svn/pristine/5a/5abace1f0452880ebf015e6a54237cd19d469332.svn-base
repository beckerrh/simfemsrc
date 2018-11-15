#ifndef __FadalightMesh_HangingNodeInfo_h
#define __FadalightMesh_HangingNodeInfo_h

#include  "FadalightMesh/geometryobject.hpp"
#include  "Alat/sparsitypattern.hpp"

/*---------------------------------------------------------*/
namespace FadalightMesh
{
  class HangingNodeInfo : public FadalightMesh::GeometryObject, public alat::SparsityPattern
  {
  public:
    ~HangingNodeInfo();
    HangingNodeInfo();
    HangingNodeInfo(const HangingNodeInfo& hangingnodeinfo);
    HangingNodeInfo& operator=(const HangingNodeInfo& hangingnodeinfo);

    std::string getClassName() const;
    int getCellNumber(int i) const;
    int getLocalSide(int i) const;
    int getNumberOfHangingNodes(int i) const;
    int getHangingNodes(int i, int in) const;
    int& getCellNumber(int i);
    int& getLocalSide(int i);
    int& getHangingNodes(int i, int in);
    void set_size(int n_hanging, int n_local_data);
    void load(std::string filename );
    void save(std::string filename, arma::file_type datatype = arma::arma_binary) const;
  };
}


#endif

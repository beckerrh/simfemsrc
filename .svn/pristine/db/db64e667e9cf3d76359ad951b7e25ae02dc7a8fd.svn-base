#ifndef __FadalightMesh_PatchInfo_h
#define __FadalightMesh_PatchInfo_h

#include  "Alat/vector.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
namespace FadalightMesh
{
class PatchInfo
{
public:
  ~PatchInfo();
  PatchInfo();
  PatchInfo( const PatchInfo& patchinfo);
  PatchInfo& operator=( const PatchInfo& patchinfo);
  std::string getClassName() const;
  PatchInfo* clone() const;

  alat::Vector<alat::armaivec> cells, edgesinner, edgesouter, nodes;

  arma::field<arma::mat> sidecoeffs;
  arma::field<arma::imat> sideindices;
};
}

/*--------------------------------------------------------------------------*/
#endif

#ifndef __FadalightMesh_CouplingSideInformation_h
#define __FadalightMesh_CouplingSideInformation_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/node.hpp"
#include  "Alat/armadillo.hpp"
#include  "Alat/vector.hpp"

/*--------------------------------------------------------------------------*/
namespace FadalightMesh
{
  class CouplingSideInformation : public alat::InterfaceBase
  {
protected:
    int _ndata;
    alat::Vector<alat::Node> _originl;
    alat::Vector<alat::Node> _originr;
    alat::armamat _hsub;
    alat::armamat _hl;
    alat::armamat _hr;
    alat::armaivec _ikl;
    alat::armaivec _ikr;
    std::string getInterfaceName() const;
    alat::armaivec _nodeidsl, _nodeidsr;

public:
    ~CouplingSideInformation();
    CouplingSideInformation();
    CouplingSideInformation( const CouplingSideInformation& couplingsideinformation);
    CouplingSideInformation& operator=( const CouplingSideInformation& couplingsideinformation);
    std::string getClassName() const;
    void set_size(int n, int dimension = 2);
    CouplingSideInformation* clone() const;
    int getNCouplingDatas() const;
    int getiKL(int icoupling) const;
    // int& getiKL(int icoupling);
    int getiKR(int icoupling) const;
    // int& getiKR(int icoupling);
    alat::Node getOriginL(int icoupling) const;
    alat::Node& getOriginL(int icoupling);
    alat::Node getOriginR(int icoupling) const;
    alat::Node& getOriginR(int icoupling);
    double getHsub(int icoupling, int idir = 0) const;
    double& getHsub(int icoupling, int idir = 0);
    double getHL(int icoupling, int idir = 0) const;
    double& getHL(int icoupling, int idir = 0);
    double getHR(int icoupling, int idir = 0) const;
    double& getHR(int icoupling, int idir = 0);
    const alat::armaivec& getiKL() const;
    alat::armaivec& getiKL();
    const alat::armaivec& getiKR() const;
    alat::armaivec& getiKR();
    const alat::armaivec& getNodeIdsL() const;
    alat::armaivec &getNodeIdsL();
    const alat::armaivec& getNodeIdsR() const;
    alat::armaivec& getNodeIdsR();
    int getNodeIdsL(int icoupling) const;
    // int& getNodeIdsL(int icoupling);
    int getNodeIdsR(int icoupling) const;
    // int& getNodeIdsR(int icoupling);
    void write(std::string filename, arma::file_type datatype = arma::arma_binary) const;
    void read(std::string filename);
  };
}

/*--------------------------------------------------------------------------*/

#endif

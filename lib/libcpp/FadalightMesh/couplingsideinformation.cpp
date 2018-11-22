#include  "FadalightMesh/couplingsideinformation.hpp"
#include <fstream>
#include  <cassert>

using namespace FadalightMesh;

/*--------------------------------------------------------------------------*/
CouplingSideInformation::~CouplingSideInformation(){}
CouplingSideInformation::CouplingSideInformation() : alat::InterfaceBase(){}
CouplingSideInformation::CouplingSideInformation( const CouplingSideInformation& couplingsideinformation) : alat::InterfaceBase(couplingsideinformation)
{
//  assert(0);
}
CouplingSideInformation& CouplingSideInformation::operator=( const CouplingSideInformation& couplingsideinformation)
{
  alat::InterfaceBase::operator=(couplingsideinformation);
  assert(0);
  return *this;
}
std::string CouplingSideInformation::getClassName() const
{
  return "CouplingSideInformation";
}
std::string CouplingSideInformation::getInterfaceName() const
{
  return "CouplingMeshInterface";
}
CouplingSideInformation* CouplingSideInformation::clone() const
{
  assert(0);
  return new CouplingSideInformation(*this);
}

/*--------------------------------------------------------------------------*/

int CouplingSideInformation::getNCouplingDatas() const
{
    return _ndata;
}

/*--------------------------------------------------------------------------*/

int CouplingSideInformation::getiKL(int icoupling) const
{
    return _ikl[icoupling];
}

/*--------------------------------------------------------------------------*/

int CouplingSideInformation::getiKR(int icoupling) const
{
   return  _ikr[icoupling];
}

/*--------------------------------------------------------------------------*/

alat::Node CouplingSideInformation::getOriginL(int icoupling) const
{
  return _originl[icoupling];
}

/*--------------------------------------------------------------------------*/

alat::Node &CouplingSideInformation::getOriginL(int icoupling)
{
  return _originl[icoupling];
}
/*--------------------------------------------------------------------------*/

alat::Node CouplingSideInformation::getOriginR(int icoupling) const
{
  return _originr[icoupling];
}

/*--------------------------------------------------------------------------*/

alat::Node &CouplingSideInformation::getOriginR(int icoupling)
{
  return _originr[icoupling];
}

/*--------------------------------------------------------------------------*/

double CouplingSideInformation::getHsub(int icoupling, int idir) const
{
  return _hsub(icoupling,idir);
}

/*--------------------------------------------------------------------------*/

double& CouplingSideInformation::getHsub(int icoupling, int idir)
{
  return _hsub(icoupling,idir);
}

/*--------------------------------------------------------------------------*/

double CouplingSideInformation::getHL(int icoupling, int idir) const
{
  return _hl(icoupling,idir);
}

/*--------------------------------------------------------------------------*/

double& CouplingSideInformation::getHL(int icoupling, int idir)
{
  return _hl(icoupling,idir);
}

/*--------------------------------------------------------------------------*/

double CouplingSideInformation::getHR(int icoupling, int idir) const
{
  return _hr(icoupling,idir);
}

/*--------------------------------------------------------------------------*/

double& CouplingSideInformation::getHR(int icoupling, int idir)
{
  return _hr(icoupling,idir);
}

/*--------------------------------------------------------------------------*/

const alat::armaivec& CouplingSideInformation::getiKL() const
{
    return _ikl;
}

/*--------------------------------------------------------------------------*/

alat::armaivec& CouplingSideInformation::getiKL()
{
    return _ikl;
}

/*--------------------------------------------------------------------------*/

const alat::armaivec& CouplingSideInformation::getiKR() const
{
    return _ikr;
}

/*--------------------------------------------------------------------------*/

alat::armaivec& CouplingSideInformation::getiKR()
{
    return _ikr;
}

/*--------------------------------------------------------------------------*/

alat::armaivec& CouplingSideInformation::getNodeIdsL()
{
    return _nodeidsl;
}

/*--------------------------------------------------------------------------*/

const alat::armaivec& CouplingSideInformation::getNodeIdsL() const
{
    return _nodeidsl;
}

/*--------------------------------------------------------------------------*/

alat::armaivec& CouplingSideInformation::getNodeIdsR()
{
    return _nodeidsr;
}

/*--------------------------------------------------------------------------*/

const alat::armaivec& CouplingSideInformation::getNodeIdsR() const
{
    return _nodeidsr;
}

/*--------------------------------------------------------------------------*/

int CouplingSideInformation::getNodeIdsL(int icoupling) const
{
    return _nodeidsl[icoupling];
}

/*--------------------------------------------------------------------------*/

int CouplingSideInformation::getNodeIdsR(int icoupling) const
{
    return _nodeidsr[icoupling];
}

/*--------------------------------------------------------------------------*/
void CouplingSideInformation::write(std::string filename, arma::file_type datatype) const
{
    std::ofstream file( filename.c_str() );
    assert( file.is_open() );
    _originl.save(file,datatype);
    _originr.save(file,datatype);
    _hsub.save(file,datatype);
    _hl.save(file,datatype);
    _hr.save(file,datatype);
    _ikl.save(file,datatype);
    _ikr.save(file,datatype);
    file.close();
}

/*--------------------------------------------------------------------------*/

void CouplingSideInformation::read(std::string filename)
{
    std::ifstream file( filename.c_str() );
    assert( file.is_open() );
    _originl.load(file);
    _originr.load(file);
    _hsub.load(file);
    _hl.load(file);
    _hr.load(file);
    _ikl.load(file);
    _ikr.load(file);
    _ndata=_ikr.size();
}
/*--------------------------------------------------------------------------*/
void CouplingSideInformation::set_size(int n, int dimension)
{
    _ndata=n;
    _originl.set_size(n);
    _originr.set_size(n);
    _hsub.set_size(n,dimension-1);
    _hl.set_size(n,dimension-1);
    _hr.set_size(n,dimension-1);
    _ikl.set_size(n);
    _ikr.set_size(n);
    _nodeidsl.set_size(n+1);
    _nodeidsr.set_size(n+1);
}

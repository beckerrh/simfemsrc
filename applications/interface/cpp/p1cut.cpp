#include  "p1cut.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Alat/matrixonevariable.hpp"
#include  "Alat/matrixallvariables.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
P1Cut::~P1Cut() {}
P1Cut::P1Cut(): P1Hansbo()
{
}
P1Cut::P1Cut( const P1Cut& P1cut): P1Hansbo(P1cut)
{
  assert(0);
}
P1Cut& P1Cut::operator=( const P1Cut& P1cut)
{
  assert(0);
  P1Hansbo::operator=(P1cut);
  return *this;
}
std::string P1Cut::getClassName() const
{
  return "P1Cut";
}

/*--------------------------------------------------------------------------*/
void P1Cut::initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp)
{
  // std::cerr << "P1Cut::initFem()\n";
  P1Hansbo::initFem(ivar, mesh, meshinfo, ncomp);
  int n = _mesh->getNNodesPerCell();
  beta.set_size(3);
  Iin.set_size(3,3);
  Iex.set_size(3,3);
  P.set_size(3,3);
  Pnc.set_size(3,3);
  Qnc.set_size(3,3);
  Iin.zeros();
  Iex.zeros();
  Id = arma::eye<arma::mat>(3,3);
}
/*--------------------------------------------------------------------------*/
void P1Cut::setCell(int iK)
{
  P1Hansbo::setCell(iK);
  int iKcut = (*_celliscut)[iK];
  if(iKcut<0)
  {
    return;
  }
  int nn = _mesh->getNNodesPerCell();
  int nnodes = _mesh->getNNodes();
  Iin.zeros(); Iex.zeros();
  for(int ii=0;ii<nn;ii++)
  {
    if((*_nodesofcutcellsisin)(ii,iKcut))
    {
      Iin(ii,ii)=1.0;
      Iex(ii,ii)=0.0;
    }
    else
    {
      Iin(ii,ii)=0.0;
      Iex(ii,ii)=1.0;
    }
  }

  // test
  // alat::armaivec indices;
  // indicesOfCell(iK, indices);
  // for(int ii=0; ii<nn;ii++)
  // {
  //   int iNin = indices[indin[ii]];
  //   int iNex = indices[index[ii]];
  //   if(iNin<_meshinfo->nnodes)
  //   {
  //     assert(iNex>=_meshinfo->nnodes);
  //     int iN = iNin;
  //     int iNcut = (*_nodeiscut)[iN];
  //     assert(iNcut+_meshinfo->nnodes==iNex);
  //     assert(not (*_cutnodesisin)[iNcut]);
  //   }
  //   else
  //   {
  //     assert(iNex<_meshinfo->nnodes);
  //     int iN = iNex;
  //     int iNcut = (*_nodeiscut)[iN];
  //     assert(iNcut+_meshinfo->nnodes==iNin);
  //     assert((*_cutnodesisin)[iNcut]);
  //   }
  // }

  // std::cerr << "indin="<<indin.t();
  // std::cerr << "index="<<index.t();
}
/*--------------------------------------------------------------------------*/
void P1Cut::computeBeta(int iK)
{
  const mesh::EdgesAndCells& edgesandcells = _mesh->getEdgesAndCells();
  int iKcut = (*_celliscut)[iK];
  assert(iKcut>=0);
  const arma::subview_col<double> normal = (*_normalsofcutcells).col(iKcut);
  int nec = _mesh->getNEdgesPerCell();
  double b=0.0;
  arma::mat B = arma::zeros<arma::mat>(2,3);
  arma::mat Bnc = arma::zeros<arma::mat>(1,3);
  int count=0;
  // Pin.zeros(); Pex.zeros();
  for(int ii=0;ii<nec;ii++)
  {
    int ie = edgesandcells._edges_of_cells(ii,iK);
    int ile = (*_edgeiscut)[ie];
    if(ile<0) continue;
    int ice = (*_cutedges)[ile];
    double cutcoeff = (*_cutcoeff)[ile];
    int i0 = edgesandcells._localnodes_of_edges_in_cells(0, ii, iK);
    int i1 = edgesandcells._localnodes_of_edges_in_cells(1, ii, iK);

    double bp = (1.0-cutcoeff);
    double ap = cutcoeff;
    double cp = ap*bp/(ap*ap+bp*bp);
    ap /= (ap*ap+bp*bp);
    bp /= (ap*ap+bp*bp);
    if(not (*_nodesofcutcellsisin)(i0,iKcut))
    {
      // assert(Pin(i0,i0)==0.0);
      // assert(Pin(i1,i1)==0.0);
      // Pin(i0,i0) += 1.0-ap;
      // Pin(i0,i1) -= cp;
      // Pin(i1,i0) -= cp;
      // Pin(i1,i1) += 1.0-bp;
    }
    else
    {
      // assert(Pex(i0,i0)==0.0);
      // assert(Pex(i1,i1)==0.0);
      // Pex(i0,i0) += 1.0-ap;
      // Pex(i0,i1) -= cp;
      // Pex(i1,i0) -= cp;
      // Pex(i1,i1) += 1.0-bp;
    }

    // conf
    B(count, i0) += (1.0-cutcoeff);
    B(count, i1) += cutcoeff;
    count++;

    // nonconf
    Bnc(0, i0) += 0.5*(1.0-cutcoeff);
    Bnc(0, i1) += 0.5*cutcoeff;

    int iN0 = _meshinfo->nodes_of_cells(i0,iK);
    int iN1 = _meshinfo->nodes_of_cells(i1,iK);

    // std::cerr << "cut " << iN0 <<  ":" << iN1 << " -> " << cutcoeff<<"\n";

    alat::armavec xc, x0, x1;
    x0 = _meshinfo->nodes.col(iN0);
    x1 = _meshinfo->nodes.col(iN1);
    xc = (1.0-cutcoeff) *x0  + cutcoeff*x1;
    b = arma::dot(normal, xc);
  }
  // std::cerr << "B="<<B<< "pB=" << arma::pinv(B);
  for(int ii=0;ii<3;ii++)
  {
    int iN = _meshinfo->nodes_of_cells(ii,iK);
    beta[ii] = arma::dot(_meshinfo->nodes.col(iN), normal) -b;
  }

  // Q = beta*beta.t()/(2.0*arma::dot(beta,beta));
  // Q = beta*beta.t()/(arma::dot(beta,beta));
  // P = arma::eye<arma::mat>(size(P)) - Q;
  // P = Id - Q;

  // std::cerr << "beta="<<beta/arma::norm(beta);
  // std::cerr << "P="<<P;
  // std::cerr << "Pin+Pex="<<Pin+Pex;
  // std::cerr << "arma::pinv(B)*B="<<arma::pinv(B)*B;
  // assert(arma::norm(P- Pin-Pex) < 1e-10);

  P = 0.5*arma::pinv(B)*B;
  Q = Id - P;
  Pnc = 0.5*arma::pinv(Bnc)*Bnc;
  Qnc = Id - Pnc;

  // std::cerr << "beta="<<beta.t();
  // int nn = _mesh->getNNodesPerCell();
  // Iin.zeros(); Iex.zeros();
  // for(int ii=0;ii<nn;ii++)
  // {
  //   if((*_nodesofcutcellsisin)(ii,iKcut))
  //   {
  //     Iin(ii,ii)=1.0;
  //     Iex(ii,ii)=0.0;
  //   }
  //   else
  //   {
  //     Iin(ii,ii)=0.0;
  //     Iex(ii,ii)=1.0;
  //   }
  // }
}

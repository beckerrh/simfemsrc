#include  "p1iim.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
P1IIM::~P1IIM() {}
P1IIM::P1IIM(): solvers::P1()
, _cutcells(NULL),_cutedges(NULL),_celliscut(NULL), _edgeiscut(NULL),_cutnodes(NULL),_cutcoeff(NULL),_nodesofcutcellsisin(NULL),_nodesofcutcells(NULL),_measuresofcutcells(NULL),_normalsofcutcells(NULL)
{
}
P1IIM::P1IIM( const P1IIM& P1cut): solvers::P1(P1cut)
{
  assert(0);
}
P1IIM& P1IIM::operator=( const P1IIM& P1cut)
{
  assert(0);
  solvers::P1::operator=(P1cut);
  return *this;
}
std::string P1IIM::getClassName() const
{
  return "P1IIM";
}

/*--------------------------------------------------------------------------*/
void P1IIM::initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp)
{
  initCutInterface(mesh);
  solvers::P1::initFem(ivar, mesh, meshinfo, ncomp);
  int n = _mesh->getNNodesPerCell();
  laplaceCut.set_size(n,n);
  femcoefin.set_size(n,n);
  femcoefex.set_size(n,n);
  Iin.set_size(n,n);
  Iex.set_size(n,n);
  Iin.zeros();
  Iex.zeros();
  Bin.set_size(n,n);
  Bex.set_size(n,n);
  C.set_size(n,n);
  P.set_size(n,n);
  laplaceP1.set_size(n,n);
}
void P1IIM::setCoefs(double kin, double kex)
{
  _kin = kin;
  _kex = kex;
}

/*--------------------------------------------------------------------------*/
void P1IIM::initCutInterface(const mesh::MeshUnitInterface* mesh)
{
  assert(mesh->geometryObjectExists(meshEnums::CutInterface));
  std::shared_ptr<const mesh::CutInterface> _cutinterface = std::dynamic_pointer_cast<const mesh::CutInterface>(mesh->getGeometryObject(meshEnums::CutInterface));
  assert(_cutinterface);
  _cutcells = &_cutinterface->getCutCells();
  _cutedges = &_cutinterface->getCutEdges();
  _celliscut = &_cutinterface->getCellIsCut();
  _edgeiscut = &_cutinterface->getEdgeIsCut();
  _cutnodes = &_cutinterface->getCutNodes();
  _cutcoeff = &_cutinterface->getCutCoeff();
  _nodesofcutcellsisin = &_cutinterface->getNodesOfCutCellsIsIn();
  _nodesofcutcells = &_cutinterface->getNodesOfCutCells();
  _measuresofcutcells = &_cutinterface->getMeasureOfCutCells();
  _normalsofcutcells = &_cutinterface->getNormalsOfCutCells();
  _cinofcutcells = &_cutinterface->getCInOfCutCells();
  _cexofcutcells = &_cutinterface->getCExOfCutCells();
  _cofinofcutcells = &_cutinterface->getCofInOfCutCells();
  _cofexofcutcells = &_cutinterface->getCofExOfCutCells();
}
/*--------------------------------------------------------------------------*/
void P1IIM::setCell(int iK)
{
  solvers::P1::setCell(iK);
  laplaceP1 = _femdata.laplace;
  iKcut = (*_celliscut)[iK];
  if(iKcut<0)
  {
    if(iKcut==-1) iKisin=true;
    else iKisin=false;
    return;
  }
  computeCutFem(iK);
}
/*--------------------------------------------------------------------------*/
void P1IIM::computeCutFem(int iK)
{
  const mesh::EdgesAndCells& edgesandcells = _mesh->getEdgesAndCells();
  iKcut = (*_celliscut)[iK];
  assert(iKcut>=0);

  int nn = _mesh->getNNodesPerCell();
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
  const arma::subview_col<double> normal = (*_normalsofcutcells).col(iKcut);
  double dnorm = arma::norm(normal);
  int nec = _mesh->getNEdgesPerCell();
  int count=0;
  Bin.zeros(); Bex.zeros(); C.zeros();
  for(int iie=0;iie<nec;iie++)
  {
    int ie = edgesandcells._edges_of_cells(iie,iK);
    int ile = (*_edgeiscut)[ie];
    if(ile<0) continue;
    double cutcoeff = (*_cutcoeff)[ile];
    int i0 = edgesandcells._localnodes_of_edges_in_cells(0, iie, iK);
    int i1 = edgesandcells._localnodes_of_edges_in_cells(1, iie, iK);
    Bin(count,i0) = (1.0-cutcoeff);
    Bin(count,i1) = cutcoeff;
    Bex(count,i0) = (1.0-cutcoeff);
    Bex(count,i1) = cutcoeff;
    count++;
  }
  assert(count==2);
  for(int ii=0;ii<3;ii++)
  {
    double dphin = arma::dot(_femdata.dphi.col(ii),normal);//dnorm;
    Bin(2, ii) = _kin*dphin;
    Bex(2, ii) = _kex*dphin;
    C(2, ii) = (_kex-_kin)*dphin;
  }

    P = arma::pinv(Bin)*Bex;
    femcoefex = arma::pinv(Iin*P + Iex).t();
    femcoefin = femcoefex*P.t();


  arma::mat A(6,6), D(6,3);
  A.submat(0,0,2,2) = Iin;
  A.submat(0,3,2,5) = Iex;
  A.submat(3,0,5,2) = Bin;
  A.submat(3,3,5,5) = -Bex;
  D.submat(0,0,2,2) = arma::zeros<arma::mat>(3,3);
  D.submat(3,0,5,2) = C;
  arma::mat X = arma::pinv(A)*D;
  // arma::mat X = arma::solve(A,D,arma::solve_opts::equilibrate	);
  femcoefin = X.submat(0,0,2,2).t() + arma::eye<arma::mat>(3,3);
  femcoefex = X.submat(3,0,5,2).t() + arma::eye<arma::mat>(3,3);

  // std::cerr << "iK=" << iK<<"\n";
  // std::cerr << "femcoefin\n" << femcoefin;
  // std::cerr << "femcoefex\n" << femcoefex << "\n";

  // for(int iie=0;iie<nec;iie++)
  // {
  //   int ie = edgesandcells._edges_of_cells(iie,iK);
  //   int ile = (*_edgeiscut)[ie];
  //   if(ile<0) continue;
  //   double cutcoeff = (*_cutcoeff)[ile];
  //   int i0 = edgesandcells._localnodes_of_edges_in_cells(0, iie, iK);
  //   int i1 = edgesandcells._localnodes_of_edges_in_cells(1, iie, iK);
  //   int iN0 = _meshinfo->nodes_of_cells(i0,iK);
  //   int iN1 = _meshinfo->nodes_of_cells(i1,iK);
  //   std::cerr << " iN=" << iN0<<" : "<< iN1<<" "<< (1.0-cutcoeff)*_meshinfo->nodes.col(iN0) + cutcoeff*_meshinfo->nodes.col(iN1);
  //   for(int ii=0;ii<nn;ii++)
  //   {
  //     double lin = (1.0-cutcoeff)*femcoefin(ii,i0) + cutcoeff*femcoefin(ii,i1);
  //     double lex = (1.0-cutcoeff)*femcoefex(ii,i0) + cutcoeff*femcoefex(ii,i1);
  //     std::cerr << "ii="<<ii<<" lr = " << (1.0-cutcoeff)*femcoefin(ii,i0) + cutcoeff*femcoefin(ii,i1)<< " => "<< (1.0-cutcoeff)*femcoefex(ii,i0) + cutcoeff*femcoefex(ii,i1)<<"\n";
  //     assert(fabs(lin-lex)<1e-10);
  //   }
  // }
  // double dnin=0.0, dnex=0.0;
  // for(int ii=0;ii<3;ii++)
  // {
  //   for(int jj=0;jj<3;jj++)
  //   {
  //     double dphin = arma::dot(_femdata.dphi.col(jj),normal)/dnorm;
  //     dnin += _kin*femcoefin(ii,jj)*dphin;
  //     dnex += _kex*femcoefex(ii,jj)*dphin;
  //   }
  //   std::cerr << "dnin="<<dnin << " dnex="<<dnex<<"\n";
  //   assert(fabs(dnin-dnex)<1e-10);
  // }
  // for(int ii=0;ii<nn;ii++)
  // {
  //   for(int jj=0;jj<nn;jj++)
  //   {
  //     double deltaij=0.0;
  //     if(jj==ii) deltaij=1.0;
  //     if((*_nodesofcutcellsisin)(jj,iKcut))
  //     {
  //       assert( fabs(femcoefin(ii,jj)- deltaij)<1e-10);
  //     }
  //     else
  //     {
  //       assert( fabs(femcoefex(ii,jj)- deltaij)<1e-10);
  //     }
  //   }
  // }

  double moc=_meshinfo->measure_of_cells[iK];
  double Kin = (*_measuresofcutcells)(0,iKcut);
  double Kex = (*_measuresofcutcells)(1,iKcut);
  double diffin = _kin * Kin/moc;
  double diffex = _kex * Kex/moc;
  laplaceCut = diffin*(femcoefin*_femdata.laplace*femcoefin.t()) + diffex*(femcoefex*_femdata.laplace*femcoefex.t());
  // std::cerr << "iK="<<iK<<" moc="<<moc << " Kin="<<Kin<<" Kex="<<Kex<< " _kin="<<_kin<<" _kex="<<_kex<<"\n";
  // std::cerr << "_femdata.laplace="<<_femdata.laplace;
  // std::cerr << "laplaceCut="<<laplaceCut;
}

/*--------------------------------------------------------------------------*/
void P1IIM::computeErrors(int iK, solvers::ErrorsMap& errormaps, const alat::armavec& uloc, const solvers::FunctionInterface& exactsolutions)
{
  setCell(iK);
  int iKcut = (*_celliscut)[iK];
  if(iKcut<0)
  {
    double k = _kin;
    if (iKcut==-2) k = _kex;
    alat::armavec uex(_ncomp), uex_x(_ncomp), uex_y(_ncomp), uex_z(_ncomp);
    arma::mat ugradex(3,_ncomp);
    int nloccell = getNPerCell(iK);
    const solvers::IntegrationFormulaInterface* IF = getFormulaErrors();
    // setCell(iK);
    for(int k = 0; k < IF->n(); k++)
    {
      const solvers::FemData& fem = referencePointWithData(IF->point(k),IF->weight(k), uloc);
      exactsolutions(uex, fem.x, fem.y, fem.z);
      assert(uex.is_finite());
      uex -= fem.u;
      if(errormaps.hasKey("L1"))
      {
        errormaps["L1"] += fem.weight*arma::norm(uex, 1);
      }
      if(errormaps.hasKey("Linf"))
      {
        for(int icomp=0;icomp<errormaps["Linf"].size();icomp++)
        {
          errormaps["Linf"][icomp] = fmax(errormaps["Linf"][icomp], fabs(uex[icomp]));
        }
      }
      if(errormaps.hasKey("L2"))
      {
        // std::cerr << fem.weight <<  " : " << uex*uex << "\n";
        errormaps["L2"] += fem.weight*(uex*uex);
      }
      if(errormaps.hasKey("H1"))
      {
        exactsolutions.x(uex_x, fem.x, fem.y, fem.z);
        exactsolutions.y(uex_y, fem.x, fem.y, fem.z);
        exactsolutions.z(uex_z, fem.x, fem.y, fem.z);
        ugradex.row(0) = uex_x;
        ugradex.row(1) = uex_y;
        ugradex.row(2) = uex_z;
        ugradex -= fem.ugrad;
        errormaps["H1"] += fem.weight*arma::dot(ugradex,ugradex);
      }
      if(errormaps.hasKey("E"))
      {
        exactsolutions.x(uex_x, fem.x, fem.y, fem.z);
        exactsolutions.y(uex_y, fem.x, fem.y, fem.z);
        exactsolutions.z(uex_z, fem.x, fem.y, fem.z);
        ugradex.row(0) = uex_x;
        ugradex.row(1) = uex_y;
        ugradex.row(2) = uex_z;
        ugradex -= fem.ugrad;
        errormaps["E"] += k*fem.weight*arma::dot(ugradex,ugradex);
      }
    }
    return;
  }
  // std::cerr << "errormaps"<<errormaps;

  const arma::subview_col<double> xin = (*_cinofcutcells).col(iKcut);
  const arma::subview_col<double> xex = (*_cexofcutcells).col(iKcut);
  const arma::subview_col<double> coeffin = (*_cofinofcutcells).col(iKcut);
  const arma::subview_col<double> coeffex = (*_cofexofcutcells).col(iKcut);

  alat::armavec uexin(_ncomp), uex_xin(_ncomp), uex_yin(_ncomp), uex_zin(_ncomp);
  alat::armavec uexex(_ncomp), uex_xex(_ncomp), uex_yex(_ncomp), uex_zex(_ncomp);
  arma::mat ugradex(3,_ncomp);
  int nloccell = getNPerCell(iK);

  alat::armaivec indices;
  indicesOfCell(iK, indices);

  // computeCutCentersOfGravity(iK);
  double Kin = (*_measuresofcutcells)(0,iKcut);
  double Kex = (*_measuresofcutcells)(1,iKcut);

  exactsolutions(uexin, xin[0], xin[1], xin[2]);
  exactsolutions.x(uex_xin, xin[0], xin[1], xin[2]);
  exactsolutions.y(uex_yin, xin[0], xin[1], xin[2]);
  exactsolutions.z(uex_zin, xin[0], xin[1], xin[2]);

  exactsolutions(uexex, xex[0], xex[1], xex[2]);
  exactsolutions.x(uex_xex, xex[0], xex[1], xex[2]);
  exactsolutions.y(uex_yex, xex[0], xex[1], xex[2]);
  exactsolutions.z(uex_zex, xex[0], xex[1], xex[2]);

  int nn = _mesh->getNNodesPerCell();
  alat::armavec uhin(_ncomp), uhex(_ncomp);
  arma::mat uhgradin(3,_ncomp), uhgradex(3,_ncomp);

  // std::cerr << "iK="<<iK << "\n";
  // std::cerr << "xin="<<xin.t();
  // std::cerr << "xex="<<xex.t();

  // alat::armavec xin2(3), xex2(3);
  // xin2.zeros(); xex2.zeros();
  // for(int ii=0;ii<nn;ii++)
  // {
  //   int iN = _meshinfo->nodes_of_cells(ii, iK);
  //   std::cerr << "iN = " << iN << " u= " << uloc(0, ii) << "\n";
  //   xin2 += coeffin[ii]*_meshinfo->nodes.col(iN);
  //   xex2 += coeffex[ii]*_meshinfo->nodes.col(iN);
  // }
  // assert( arma::norm(xin-xin2) < 1e-10);
  // assert( arma::norm(xex-xex2) < 1e-10);

  uhin.zeros(); uhgradin.zeros();
  uhex.zeros(); uhgradex.zeros();
  alat::armavec lambdahatin = femcoefin*coeffin;
  alat::armavec lambdahatex = femcoefex*coeffex;
  // arma::mat dlambdahatin, dlambdahatex;
  arma::mat dlambdahatin = _femdata.dphi*femcoefin.t();
  arma::mat dlambdahatex = _femdata.dphi*femcoefex.t();
  for(int ii=0;ii<nn;ii++)
  {
    // std::cerr << "\t lambdai=" << lambdahatin[ii] << " coeffin=" << coeffin[ii] << " ui="<< uloc(0, ii) <<"\n";
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      uhin[icomp] += lambdahatin[ii]*uloc(icomp, ii);
      uhex[icomp] += lambdahatex[ii]*uloc(icomp, ii);
      uhgradin.col(icomp) += dlambdahatin.col(ii)*uloc(icomp, ii);
      uhgradex.col(icomp) += dlambdahatex.col(ii)*uloc(icomp, ii);
    }
  }
  // std::cerr << "uexin= "<<uexin[0] << " uexex= "<<uexex[0];
  // std::cerr << " uhin= "<<uhin[0]<< " uhex= "<<uhex[0] << "\n";

  // std::cerr << "uex_xin= "<<uex_xin[0] << " uex_yin= "<<uex_yin[0];
  // std::cerr << " == uhgradin= "<< alat::armavec(uhgradin.col(0)).t();
  //
  // std::cerr << "uex_xex= "<<uex_xex[0] << " uex_uex= "<<uex_yex[0];
  // std::cerr << " == uhgradex= "<< alat::armavec(uhgradex.col(0)).t();

  uexin -= uhin;
  uexex -= uhex;
  if(errormaps.hasKey("Linf"))
  {
    for(int icomp=0;icomp<errormaps["Linf"].size();icomp++)
    {
      errormaps["Linf"][icomp] = fmax(errormaps["Linf"][icomp], fabs(uexin[icomp]));
      errormaps["Linf"][icomp] = fmax(errormaps["Linf"][icomp], fabs(uexex[icomp]));
    }
  }
  if(errormaps.hasKey("L1"))
  {
    errormaps["L1"] += Kin*arma::norm(uexin, 1);
    errormaps["L1"] += Kex*arma::norm(uexex, 1);
  }
  if(errormaps.hasKey("L2"))
  {
    // std::cerr << fem.weight <<  " : " << uex*uex << "\n";
    errormaps["L2"] += Kin*(uexin*uexin);
    errormaps["L2"] += Kex*(uexex*uexex);
  }
  if(errormaps.hasKey("H1"))
  {
    ugradex.row(0) = uex_xin;
    ugradex.row(1) = uex_yin;
    ugradex.row(2) = uex_zin;
    ugradex -= uhgradin;
    errormaps["H1"] += Kin*arma::dot(ugradex,ugradex);
    ugradex.row(0) = uex_xex;
    ugradex.row(1) = uex_yex;
    ugradex.row(2) = uex_zex;
    ugradex -= uhgradex;
    errormaps["H1"] += Kex*arma::dot(ugradex,ugradex);
  }
  if(errormaps.hasKey("E"))
  {
    ugradex.row(0) = uex_xin;
    ugradex.row(1) = uex_yin;
    ugradex.row(2) = uex_zin;
    ugradex -= uhgradin;
    errormaps["E"] += _kin*Kin*arma::dot(ugradex,ugradex);
    ugradex.row(0) = uex_xex;
    ugradex.row(1) = uex_yex;
    ugradex.row(2) = uex_zex;
    ugradex -= uhgradex;
    errormaps["E"] += _kex*Kex*arma::dot(ugradex,ugradex);
  }
}

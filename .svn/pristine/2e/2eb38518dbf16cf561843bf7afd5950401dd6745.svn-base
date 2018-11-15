#include  "cr1iim.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Alat/matrixonevariable.hpp"
#include  "Alat/matrixallvariables.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
CR1IIM::~CR1IIM() {}
CR1IIM::CR1IIM(): solvers::P1(), _cutcells(NULL),_cutedges(NULL),_celliscut(NULL), _edgeiscut(NULL),_cutnodes(NULL),_cutcoeff(NULL),_nodesofcutcellsisin(NULL),_nodesofcutcells(NULL),_measuresofcutcells(NULL),_normalsofcutcells(NULL)
{
}
CR1IIM::CR1IIM( const CR1IIM& P1cut): solvers::P1(P1cut)
{
  assert(0);
}
CR1IIM& CR1IIM::operator=( const CR1IIM& P1cut)
{
  assert(0);
  solvers::P1::operator=(P1cut);
  return *this;
}
std::string CR1IIM::getClassName() const
{
  return "CR1IIM";
}

/*--------------------------------------------------------------------------*/
void CR1IIM::initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp)
{
  initCutInterface(mesh);
  solvers::P1::initFem(ivar, mesh, meshinfo, ncomp);
  int n = _mesh->getNNodesPerCell();
  laplaceCut.set_size(n,n);
  femcoefin.set_size(n,n);
  femcoefex.set_size(n,n);
  Iin.set_size(n,n);
  Iex.set_size(n,n);
  Bin.set_size(n,n);
  Bex.set_size(n,n);
  P.set_size(n,n);
  laplaceP1.set_size(n,n);
  phi_P1.set_size(n);
  dphi_P1.set_size(3,n);
}
void CR1IIM::setCoefs(double kin, double kex)
{
  _kin = kin;
  _kex = kex;
}

/*--------------------------------------------------------------------------*/
void CR1IIM::initCutInterface(const mesh::MeshUnitInterface* mesh)
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
int CR1IIM::getN() const
{
  return _mesh->getNSides();
}
int CR1IIM::getNPerCell(int iK) const
{
  return _mesh->getNSidesPerCell();
}
void CR1IIM::indicesOfCell(int iK, alat::armaivec& indices) const
{
  int n = _mesh->getNSidesPerCell();
  indices.set_size(n);
  for(int ii=0;ii<n;ii++)
  {
    indices[ii] = _meshinfo->sides_of_cells(ii, iK);
  }
}
/*--------------------------------------------------------------------------*/
const solvers::FemData& CR1IIM::referencePoint(const alat::Node& vhat, double weight)
{
  _femdata.phi.zeros();
  P1::referencePoint(vhat, weight);
  phi_P1 = _femdata.phi;
  _femdata.phi.zeros();
  for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
  {
    for(int jj=0;jj<_meshinfo->nnodespercell;jj++)
    {
      _femdata.phi[ii] += femcoefin(ii,jj)*phi_P1[jj];
    }
  }
  return _femdata;
}

/*--------------------------------------------------------------------------*/
void CR1IIM::setCell(int iK)
{
  solvers::P1::setCell(iK);
  dphi_P1 = _femdata.dphi;
  laplaceP1 = _femdata.laplace;
  iKcut = (*_celliscut)[iK];
  if(iKcut<0)
  {
    if(iKcut==-1) iKisin=true;
    else iKisin=false;
    const mesh::EdgesAndCells& edgesandcells = _mesh->getEdgesAndCells();
    Iin.zeros();
    int nec = _mesh->getNEdgesPerCell();
    for(int iie=0;iie<nec;iie++)
    {
      int ie = edgesandcells._edges_of_cells(iie,iK);
      int i0 = edgesandcells._localnodes_of_edges_in_cells(0, iie, iK);
      int i1 = edgesandcells._localnodes_of_edges_in_cells(1, iie, iK);
      Iin(i0,iie) = 0.5;
      Iin(i1,iie) = 0.5;
    }
    femcoefin = arma::inv(Iin);
    _femdata.dphi.zeros();
    for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
    {
      for(int jj=0;jj<_meshinfo->nnodespercell;jj++)
      {
        _femdata.phi[ii] += femcoefin(ii,jj)*phi_P1[jj];
        _femdata.dphi.col(ii) += femcoefin(ii,jj)*dphi_P1.col(jj);
      }
    }
    // std::cerr << "dphi_P1="<<dphi_P1;
    // std::cerr << "_femdata.dphi="<<_femdata.dphi;
  }
  _femdata.laplace = femcoefin*laplaceP1*femcoefin.t();
  if(iKcut<0)
  {
    return;
  }
  computeCutFem(iK);
}

/*--------------------------------------------------------------------------*/
void CR1IIM::computeCutFem(int iK)
{
  const mesh::EdgesAndCells& edgesandcells = _mesh->getEdgesAndCells();
  iKcut = (*_celliscut)[iK];

  const arma::subview_col<double> normal = (*_normalsofcutcells).col(iKcut);
  double dnorm = arma::norm(normal);
  int count=0;
  Iin.zeros(); Iex.zeros();
  Bin.zeros(); Bex.zeros();
  int nec = _mesh->getNEdgesPerCell();
  for(int iie=0;iie<nec;iie++)
  {
    int ie = edgesandcells._edges_of_cells(iie,iK);
    int ile = (*_edgeiscut)[ie];
    int i0 = edgesandcells._localnodes_of_edges_in_cells(0, iie, iK);
    int i1 = edgesandcells._localnodes_of_edges_in_cells(1, iie, iK);
    if(ile<0)
    {
      if((*_nodesofcutcellsisin)(i0,iKcut))
      {
        assert((*_nodesofcutcellsisin)(i1,iKcut));
        // Iin(i0,iie) = 0.5;
        // Iin(i1,iie) = 0.5;
        Iin(iie,i0) = 0.5;
        Iin(iie,i1) = 0.5;
      }
      else
      {
        assert(not (*_nodesofcutcellsisin)(i1,iKcut));
        // Iex(i0,iie) = 0.5;
        // Iex(i1,iie) = 0.5;
        Iex(iie,i0) = 0.5;
        Iex(iie,i1) = 0.5;
      }
      continue;
    }
    double cutcoeff = (*_cutcoeff)[ile];
    Bin(count,i0) = (1.0-cutcoeff);
    Bin(count,i1) = cutcoeff;
    Bex(count,i0) = (1.0-cutcoeff);
    Bex(count,i1) = cutcoeff;
    count++;
    if((*_nodesofcutcellsisin)(i0,iKcut))
    {
      // Iin(i0,iie) = cutcoeff*(1.0-0.5*cutcoeff);
      // Iin(i1,iie) = 0.5*cutcoeff*cutcoeff;
      // Iex(i0,iie) = 0.5*(1.0-cutcoeff)*(1.0-cutcoeff);
      // Iex(i1,iie) = 0.5*(1.0-cutcoeff)*(1.0+cutcoeff);
      Iin(iie, i0) = cutcoeff*(1.0-0.5*cutcoeff);
      Iin(iie, i1) = 0.5*cutcoeff*cutcoeff;
      Iex(iie, i0) = 0.5*(1.0-cutcoeff)*(1.0-cutcoeff);
      Iex(iie, i1) = 0.5*(1.0-cutcoeff)*(1.0+cutcoeff);
    }
    else
    {
      // Iex(i0,iie) = cutcoeff*(1.0-0.5*cutcoeff);
      // Iex(i1,iie) = 0.5*cutcoeff*cutcoeff;
      // Iin(i0,iie) = 0.5*(1.0-cutcoeff)*(1.0-cutcoeff);
      // Iin(i1,iie) = 0.5*(1.0-cutcoeff)*(1.0+cutcoeff);
      Iex(iie, i0) = cutcoeff*(1.0-0.5*cutcoeff);
      Iex(iie, i1) = 0.5*cutcoeff*cutcoeff;
      Iin(iie, i0) = 0.5*(1.0-cutcoeff)*(1.0-cutcoeff);
      Iin(iie, i1) = 0.5*(1.0-cutcoeff)*(1.0+cutcoeff);
    }
  }
  assert(count==2);
  for(int ii=0;ii<3;ii++)
  {
    double dphin = arma::dot(_femdata.dphi.col(ii),normal)/dnorm;
    Bin(2, ii) = _kin*dphin;
    Bex(2, ii) = _kex*dphin;
  }

  // arma::mat A(6,6), D(6,3);
  // A.submat(0,0,2,2) = Iin.t();
  // A.submat(0,3,2,5) = Iex.t();
  // A.submat(3,0,5,2) = Bin;
  // A.submat(3,3,5,5) = -Bex;
  // D.submat(0,0,2,2) = arma::eye<arma::mat>(3,3);
  // D.submat(3,0,5,2) = arma::zeros<arma::mat>(3,3);
  // // arma::mat femcoefin2(3,3), femcoefex2(3,3);
  // arma::mat X = arma::inv(A)*D;
  // femcoefin = X.submat(0,0,2,2).t();
  // femcoefex = X.submat(3,0,5,2).t();

  P = arma::pinv(Bin)*Bex;
  femcoefex = arma::pinv(Iin*P + Iex).t();
  femcoefin = femcoefex*P.t();

  double moc=_meshinfo->measure_of_cells[iK];
  double Kin = (*_measuresofcutcells)(0,iKcut);
  double Kex = (*_measuresofcutcells)(1,iKcut);
  double diffin = _kin * Kin/moc;
  double diffex = _kex * Kex/moc;
  laplaceCut = diffin*(femcoefin*laplaceP1*femcoefin.t()) + diffex*(femcoefex*laplaceP1*femcoefex.t());
}
/*--------------------------------------------------------------------------*/
void CR1IIM::computeErrors(int iK, solvers::ErrorsMap& errormaps, const arma::mat& uloc, const solvers::FunctionInterface& exactsolutions)
{
  setCell(iK);
  int iKcut = (*_celliscut)[iK];
  if(iKcut<0)
  {
    double k = _kin;
    if ((*_celliscut)[iK]==-2) k = _kex;

    arma::vec uex(_ncomp), uex_x(_ncomp), uex_y(_ncomp), uex_z(_ncomp);
    arma::mat ugradex(3,_ncomp);
    int nloccell = getNPerCell(iK);
    const solvers::IntegrationFormulaInterface* IF = getFormulaErrors();
    setCell(iK);
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
  const arma::subview_col<double> xin = (*_cinofcutcells).col(iKcut);
  const arma::subview_col<double> xex = (*_cexofcutcells).col(iKcut);
  const arma::subview_col<double> coeffin = (*_cofinofcutcells).col(iKcut);
  const arma::subview_col<double> coeffex = (*_cofexofcutcells).col(iKcut);

  arma::vec uexin(_ncomp), uex_xin(_ncomp), uex_yin(_ncomp), uex_zin(_ncomp);
  arma::vec uexex(_ncomp), uex_xex(_ncomp), uex_yex(_ncomp), uex_zex(_ncomp);
  arma::mat ugradex(3,_ncomp);
  int nloccell = getNPerCell(iK);

  alat::armaivec indices;
  indicesOfCell(iK, indices);

  computeCutFem(iK);
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
  arma::vec uhin(_ncomp), uhex(_ncomp);
  arma::mat uhgradin(3,_ncomp), uhgradex(3,_ncomp);


  uhin.zeros(); uhgradin.zeros();
  uhex.zeros(); uhgradex.zeros();
  arma::vec lambdahatin = femcoefin*coeffin;
  arma::vec lambdahatex = femcoefex*coeffex;
  // arma::mat dlambdahatin, dlambdahatex;
  arma::mat dlambdahatin = dphi_P1*femcoefin.t();
  arma::mat dlambdahatex = dphi_P1*femcoefex.t();
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

  // std::cerr << "xin="<<xin.t();
  // std::cerr << "xex="<<xex.t();
  // std::cerr << "uexin= "<<uexin[0] << " uexex= "<<uexex[0];
  // std::cerr << " uhin= "<<uhin[0]<< " uhex= "<<uhex[0] << "\n";

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
    errormaps["H1"] += _kex*Kex*arma::dot(ugradex,ugradex);
  }
}
/*--------------------------------------------------------------------------*/
void CR1IIM::strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const
{
  for(alat::IntSet::const_iterator p= dircolors.begin(); p!=dircolors.end();p++)
  {
    int color = *p;
    const alat::armaimat& cells_on_bdry = _meshinfo->bdrymesheunitsmap[color].getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      int iK = cells_on_bdry(0,i);
      int iS = cells_on_bdry(1,i);
      int iil = cells_on_bdry(2,i);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        int index = icomp*_meshinfo->nsides + iS;
        A.get(ivar,ivar)->rowIdentity(index);
      }
    }
  }
}
void CR1IIM::strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const
{
  alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
  for(alat::IntSet::const_iterator p= dircolors.begin(); p!=dircolors.end();p++)
  {
    int color = *p;
    const alat::armaimat& cells_on_bdry = _meshinfo->bdrymesheunitsmap[color].getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      int iK = cells_on_bdry(0,i);
      int iS = cells_on_bdry(1,i);
      int iil = cells_on_bdry(2,i);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        int index = icomp*_meshinfo->nsides + iS;
        (*uv)[index] = 0.0;
      }
    }
  }
}
void CR1IIM::strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const
{
  alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
  arma::vec udir(_ncomp);
  for(alat::IntSet::const_iterator p= dircolors.begin(); p!=dircolors.end();p++)
  {
    int color = *p;
    const alat::armaimat& cells_on_bdry = _meshinfo->bdrymesheunitsmap[color].getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      int iK = cells_on_bdry(0,i);
      int iS = cells_on_bdry(1,i);
      int iil = cells_on_bdry(2,i);
      alat::Node xS = _mesh->getNodeOfSide(iS);
      dirichlet(udir, xS.x(), xS.y(), xS.z());
      // std::cerr << "iS="<<iS << " xy="<<xS.x()<<" "<<xS.y()<<" "<<udir.t();
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        int index = icomp*_meshinfo->nsides + iS;
        (*uv)[index] = udir[icomp];
      }
    }
  }
}

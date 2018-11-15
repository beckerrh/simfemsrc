#include  "p1hansbo.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Alat/matrixonevariable.hpp"
#include  "Alat/matrixallvariables.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
P1Hansbo::~P1Hansbo() {}
P1Hansbo::P1Hansbo(): solvers::P1(), _cutcells(NULL),_cutedges(NULL),_celliscut(NULL), _edgeiscut(NULL),_cutnodes(NULL),_cutcoeff(NULL),_nodesofcutcellsisin(NULL),_nodesofcutcells(NULL),_measuresofcutcells(NULL),_normalsofcutcells(NULL)
{
}
P1Hansbo::P1Hansbo( const P1Hansbo& P1cut): solvers::P1(P1cut)
{
  assert(0);
}
P1Hansbo& P1Hansbo::operator=( const P1Hansbo& P1cut)
{
  assert(0);
  solvers::P1::operator=(P1cut);
  return *this;
}
std::string P1Hansbo::getClassName() const
{
  return "P1Hansbo";
}

/*--------------------------------------------------------------------------*/
void P1Hansbo::initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp)
{
  // std::cerr << "P1Hansbo::initFem()\n";
  solvers::P1::initFem(ivar, mesh, meshinfo, ncomp);
  initCutInterface(mesh);
  int n = _mesh->getNNodesPerCell();
  indin.set_size(n);
  index.set_size(n);
}
void P1Hansbo::setCoefs(double kin, double kex)
{
  _kin = kin;
  _kex = kex;
}

/*--------------------------------------------------------------------------*/
void P1Hansbo::initCutInterface(const mesh::MeshUnitInterface* mesh)
{
  // std::cerr << "P1Hansbo::initCutInterface()\n";
  assert(mesh->geometryObjectExists(meshEnums::CutInterface));
  std::shared_ptr<const mesh::CutInterface> _cutinterface = std::dynamic_pointer_cast<const mesh::CutInterface>(mesh->getGeometryObject(meshEnums::CutInterface));
  assert(_cutinterface);
  _cutcells = &_cutinterface->getCutCells();
  _cutedges = &_cutinterface->getCutEdges();
  _celliscut = &_cutinterface->getCellIsCut();
  _edgeiscut = &_cutinterface->getEdgeIsCut();
  _nodeiscut = &_cutinterface->getNodeIsCut();
  _cutnodesisin = &_cutinterface->getCutNodeIsIn();
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
int P1Hansbo::getN() const
{
  return _mesh->getNNodes() + _cutnodes->size();
}
/*--------------------------------------------------------------------------*/
int P1Hansbo::getNPerCell(int iK) const
{
  int n = _mesh->getNNodesPerCell();
  if(iK==-1) return n;
  if((*_celliscut)[iK]>=0) {return 2*n;}
  return n;
}
/*--------------------------------------------------------------------------*/
void P1Hansbo::indicesOfCell(int iK, alat::armaivec& indices) const
{
  int n = _mesh->getNNodesPerCell();
  int iKcut = (*_celliscut)[iK];
  if(iKcut>=0)
  {
    int nnodes = _mesh->getNNodes();
    indices.set_size(2*n);
    for(int ii=0;ii<n;ii++)
    {
      int iN = _meshinfo->nodes_of_cells(ii, iK);
      int iNcut = (*_nodeiscut)[iN];

      // indices[ii] = iN;
      // indices[n+ii] = nnodes + iNcut;

      if((*_cutnodesisin)[iNcut])
      {
        assert((*_nodesofcutcellsisin)(ii,iKcut));
        indices[ii] = iN;
        indices[n+ii] = nnodes + iNcut;
      }
      else
      {
        indices[ii] = nnodes + iNcut;
        indices[n+ii] = iN;
      }

    }
  }
  else
  {
    indices.set_size(n);
    for(int ii=0;ii<n;ii++)
    {
      indices[ii] = _meshinfo->nodes_of_cells(ii, iK);
    }
  }
}
/*--------------------------------------------------------------------------*/
void P1Hansbo::setCell(int iK)
{
  solvers::P1::setCell(iK);
  int iKcut = (*_celliscut)[iK];
  if(iKcut<0)
  {
    // if(iKcut==-1) iKisin=true;
    // else iKisin=false;
    return;
  }
  int nn = _mesh->getNNodesPerCell();
  int nnodes = _mesh->getNNodes();
  for(int ii=0;ii<nn;ii++)
  {
    indin[ii] = ii;
    index[ii] = nn + ii;
  }
}
/*--------------------------------------------------------------------------*/
void P1Hansbo::toP1(alat::VectorOneVariableInterface* uP1, const alat::VectorOneVariableInterface* u)
{
  alat::VectorOneVariable* uP1v = dynamic_cast<alat::VectorOneVariable*>(uP1); assert(uP1v);
  const alat::VectorOneVariable* uv = dynamic_cast<const alat::VectorOneVariable*>(u); assert(uv);
  // int nnodes = _meshinfo->nnodes;
  uP1v->copyFrom(*uv);
  // solvers::P1::toP1(uP1, u);
}
/*--------------------------------------------------------------------------*/
void P1Hansbo::fromP1(alat::VectorOneVariableInterface* u, const alat::VectorOneVariableInterface* uP1)
{
  assert(0);
  solvers::P1::fromP1(u, uP1);
}
/*--------------------------------------------------------------------------*/
void P1Hansbo::computeErrors(int iK, solvers::ErrorsMap& errormaps, const arma::mat& uloc, const solvers::FunctionInterface& exactsolutions)
{
  setCell(iK);
  int iKcut = (*_celliscut)[iK];
  int nloccell = getNPerCell(iK);
  if(iKcut<0)
  {
    double k = _kin;
    if (iKcut==-2) k = _kex;
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

  alat::armaivec indices;
  indicesOfCell(iK, indices);

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

  assert(uloc.n_cols==2*nn);
  uhin.zeros(), uhgradin.zeros();
  uhex.zeros(), uhgradex.zeros();
  for(int ii=0;ii<nn;ii++)
  {
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      uhin[icomp] += coeffin[ii]*uloc(icomp, indin[ii]);
      uhex[icomp] += coeffex[ii]*uloc(icomp, index[ii]);
      uhgradin.col(icomp) += _femdata.dphi.col(ii)*uloc(icomp, indin[ii]);
      uhgradex.col(icomp) += _femdata.dphi.col(ii)*uloc(icomp, index[ii]);
    }
  }
  // std::cerr  << "coeffin= " << coeffin.t()<< " coeffex= " << coeffex.t() << "\n";
  // std::cerr  << "uexin= " << uexin<< " uhin= " << uhin;
  // std::cerr  << "uexex= " << uexex<< " uhex= " << uhex;
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
/*--------------------------------------------------------------------------*/
alat::armaivec P1Hansbo::getBdryNodes(int i, const alat::armaimat& cells_on_bdry)const
{
  int iK = cells_on_bdry(0,i);
  int iS = cells_on_bdry(1,i);
  int iil = cells_on_bdry(2,i);
  int iKcut = (*_celliscut)[iK];
  int bdrynodessize = _meshinfo->nnodesperside;
  // if(iKcut<0)
  {
    alat::armaivec bdrynodes(bdrynodessize);
    for(int ii=0;ii<bdrynodessize;ii++)
    {
      bdrynodes[ii] = _meshinfo->nodes_of_sides(ii,iS);
    }
    return bdrynodes;
  }
  alat::armaivec bdrynodes(2*bdrynodessize);
  int nnodes = _meshinfo->nnodes;
  for(int ii=0;ii<bdrynodessize;ii++)
  {
    int iN = _meshinfo->nodes_of_sides(ii,iS);
    bdrynodes[ii] = iN;
    bdrynodes[bdrynodessize+ii] = nnodes + (*_nodeiscut)[iN];
  }
  return bdrynodes;
  if(iKcut>=0)
  {
    int nnodes = _meshinfo->nnodes;
    for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
    {
      for(int jj=0;jj<_meshinfo->nnodesperside;jj++)
      {
        if(bdrynodes[jj] == _meshinfo->nodes_of_cells(ii,iK))
        {
          bdrynodes[_meshinfo->nnodesperside+jj] = nnodes + (*_nodesofcutcells)(ii,iKcut);
        }
      }
    }
    // std::cerr << "iK="<<iK << " bdrynodes="<<bdrynodes.t();
  }
  return bdrynodes;
}

/*--------------------------------------------------------------------------*/
void P1Hansbo::strongDirichlet(int ivar, alat::MatrixAllVariables& A, const alat::IntSet& dircolors)const
{
  int N = getN();
  for(alat::IntSet::const_iterator p= dircolors.begin(); p!=dircolors.end();p++)
  {
    int color = *p;
    const alat::armaimat& cells_on_bdry = _meshinfo->bdrymesheunitsmap[color].getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      alat::armaivec bdrynodes = getBdryNodes(i, cells_on_bdry);
      for(int ii=0;ii<bdrynodes.size();ii++)
      {
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          int index = icomp*N + bdrynodes[ii];
          for(int jvar=0;jvar<A.m();jvar++)
          {
            if(jvar==ivar)
            {
              A.get(ivar,jvar)->rowIdentity(index);
            }
            else
            {
              A.get(ivar,jvar)->rowZero(index);
            }
          }
        }
      }
    }
  }
}
/*--------------------------------------------------------------------------*/
void P1Hansbo::strongDirichletZero(alat::VectorOneVariableInterface* u, const alat::IntSet& dircolors)const
{
  alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
  int N = getN();
  for(alat::IntSet::const_iterator p= dircolors.begin(); p!=dircolors.end();p++)
  {
    int color = *p;
    const alat::armaimat& cells_on_bdry = _meshinfo->bdrymesheunitsmap[color].getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      alat::armaivec bdrynodes = getBdryNodes(i, cells_on_bdry);
      for(int ii=0;ii<bdrynodes.size();ii++)
      {
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          int index = icomp*N + bdrynodes[ii];
          (*uv)[index] = 0.0;
        }
      }
    }
  }
}
/*--------------------------------------------------------------------------*/
void P1Hansbo::strongDirichlet(alat::VectorOneVariableInterface* u, const solvers::DirichletInterface& dirichlet, const alat::IntSet& dircolors)const
{
  alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
  int N = getN();
  arma::vec udir(_ncomp);
  for(alat::IntSet::const_iterator p= dircolors.begin(); p!=dircolors.end();p++)
  {
    int color = *p;
    const alat::armaimat& cells_on_bdry = _meshinfo->bdrymesheunitsmap[color].getCellsOnBdryOfPlain();
    for(int i = 0; i < cells_on_bdry.n_cols; i++)
    {
      int iS = cells_on_bdry(1,i);
      alat::armaivec bdrynodes = getBdryNodes(i, cells_on_bdry);
      for(int ii=0;ii<_meshinfo->nnodesperside;ii++)
      {
        int iN = _meshinfo->nodes_of_sides(ii,iS);
        dirichlet(udir, _meshinfo->nodes(0,iN), _meshinfo->nodes(1,iN), _meshinfo->nodes(2,iN));
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          int index = icomp*N + bdrynodes[ii];
          (*uv)[index] = udir[icomp];
        }
        if(bdrynodes.size()==2*_meshinfo->nnodesperside)
        {
          for(int icomp=0;icomp<_ncomp;icomp++)
          {
            int index = icomp*N + bdrynodes[_meshinfo->nnodesperside+ii];
            (*uv)[index] = udir[icomp];
          }
        }
      }
    }
  }
}

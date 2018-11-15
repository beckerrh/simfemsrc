#include  "pdepartcutinterface.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Alat/vectorallvariables.hpp"
#include  "Alat/matrixallvariables.hpp"

/*--------------------------------------------------------------------------*/
PdePartCutInterface::~PdePartCutInterface() {}
PdePartCutInterface::PdePartCutInterface(alat::StringList vars): solvers::PdePartWithIntegration(vars){}
PdePartCutInterface::PdePartCutInterface( const PdePartCutInterface& pdepartwithfemtraditional): solvers::PdePartWithIntegration(pdepartwithfemtraditional)
{
  assert(0);
}
PdePartCutInterface& PdePartCutInterface::operator=( const PdePartCutInterface& pdepartwithfemtraditional)
{
  assert(0);
  solvers::PdePartWithIntegration::operator=(pdepartwithfemtraditional);
  return *this;
}
std::string PdePartCutInterface::getClassName() const
{
  return "PdePartCutInterface";
}
/*--------------------------------------------------------------------------*/
const solver_options::pdepart::opts PdePartCutInterface::setOptions()
{
  return solver_options::pdepart::cell;
}
/*--------------------------------------------------------------------------*/
void PdePartCutInterface::setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters)
{
  solvers::PdePartWithIntegration::setData(var2index, parameters);

  _ivar = 0;
  _ncomp = (*_fems)[_ivar]->getNcomp();
  _nlocal = (*_fems)[_ivar]->getNPerCell();
  _localmodel = dynamic_cast<const Model*>(_model);
  assert(_localmodel);
  _kin = _localmodel->_kin;
  _kex = _localmodel->_kex;

  assert(_mesh->geometryObjectExists(meshEnums::CutInterface));
  std::shared_ptr<const mesh::CutInterface> _cutinterface = std::dynamic_pointer_cast<const mesh::CutInterface>(_mesh->getGeometryObject(meshEnums::CutInterface));
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


  indicesofcutedge.resize(4,_cutedges->size());
  coefsofcutedge.resize(4,_cutedges->size());
  coefsofcutedge.zeros();
  for(int ile=0;ile<_cutedges->size();ile++)
  {
    int ie = (*_cutedges)[ile];
    double cutcoeff = (*_cutcoeff)[ile];
    double d0 = (1.0-cutcoeff);
    double d1 = cutcoeff;
    int i0 = _meshinfo->nodes_of_edges(0,ie);
    int i1 = _meshinfo->nodes_of_edges(1,ie);
    // std::cerr << "i0 i1 " << i0 << " " << i1 << " -> " << d0 << " " << d1 << "\n";
    int ic0 = (*_nodeiscut)[i0];
    int ic1 = (*_nodeiscut)[i1];
    assert(ic0>=0);
    assert(ic1>=0);
    // coefsofcutedge(0, ile) = d0;
    // coefsofcutedge(1, ile) = d1;
    // coefsofcutedge(2, ile) = -d0;
    // coefsofcutedge(3, ile) = -d1;
    if((*_cutnodesisin)[ic0])
    {
      indicesofcutedge(0, ile) = _meshinfo->nnodes+ic0;
      indicesofcutedge(1, ile) = i1;
      indicesofcutedge(2, ile) = i0;
      indicesofcutedge(3, ile) = _meshinfo->nnodes+ic1;
    }
    else
    {
      indicesofcutedge(0, ile) = i0;
      indicesofcutedge(1, ile) = _meshinfo->nnodes+ic1;
      indicesofcutedge(2, ile) = _meshinfo->nnodes+ic0;
      indicesofcutedge(3, ile) = i1;
    }
  }

  for(int iKcut=0;iKcut<_cutcells->size();iKcut++)
  {
    int iK = (*_cutcells)[iKcut];
    const arma::subview_col<double> normal = (*_normalsofcutcells).col(iKcut);
    double gammalength = arma::norm(normal);
    int nec = _mesh->getNEdgesPerCell();
    // assert(nec==3);
    const mesh::EdgesAndCells& edgesandcells = _mesh->getEdgesAndCells();
    for(int ii=0;ii<nec;ii++)
    {
      int ie = edgesandcells._edges_of_cells(ii,iK);
      int ile = (*_edgeiscut)[ie];
      if(ile<0) continue;
      double cutcoeff = (*_cutcoeff)[ile];

      double d0 = 0.5*gammalength*(1.0-cutcoeff);
      double d1 = 0.5*gammalength*cutcoeff;
      coefsofcutedge(0, ile) += d0;
      coefsofcutedge(1, ile) += d1;
      coefsofcutedge(2, ile) += -d0;
      coefsofcutedge(3, ile) += -d1;
    }
  }
}


/*--------------------------------------------------------------------------*/
void PdePartCutInterface::computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  bool celliscut = ((*_celliscut)[iK]>=0);
  if(celliscut)
  {
    computeRhsCellCut(iK, floc, uloc);
    return;
  }

  arma::vec f(_ncomp);
  double moc=_meshinfo->measure_of_cells[iK];
  alat::Node xK = _mesh->getNodeOfCell(iK);
  _application->getRightHandSide(_ivar)(f, xK.x(), xK.y(), xK.z());
  double scalemass = moc/(_meshinfo->dim+1.0);
  for(int ii=0; ii<_nlocal;ii++)
  {
    for(int icomp=0;icomp<_ncomp;icomp++)
    {
      floc[_ivar](icomp,ii) += f[icomp]*scalemass;
    }
  }
}

/*--------------------------------------------------------------------------*/
void PdePartCutInterface::computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc) const
{
  bool celliscut = ((*_celliscut)[iK]>=0);
  if(celliscut)
  {
    computeResidualCellCut(iK, floc, uloc);
    return;
  }
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  bool cellisin = (*_celliscut)[iK]==-1;
  double diff = _localmodel->_kin;
  if(not cellisin) {diff=_localmodel->_kex;}
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iS = _meshinfo->sides_of_cells(ii,iK);
    for(int jj=0; jj<_nlocal;jj++)
    {
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        floc[_ivar](icomp,ii) += diff*fem.laplace(ii,jj)*uloc[_ivar](icomp, jj);
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void PdePartCutInterface::computeMatrixCell(int iK, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const
{
  bool celliscut = ((*_celliscut)[iK]>=0);
  if(celliscut)
  {
    computeMatrixCellCut(iK, mat, mat_i, mat_j, uloc);
    return;
  }
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  bool cellisin = (*_celliscut)[iK]==-1;
  double diff = _localmodel->_kin;
  if(not cellisin) {diff=_localmodel->_kex;}
  int count=0;
  for(int ii=0; ii<_nlocal;ii++)
  {
    int iS = _meshinfo->sides_of_cells(ii,iK);
    for(int jj=0; jj<_nlocal;jj++)
    {
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        for(int jcomp=0;jcomp<_ncomp;jcomp++)
        {
          if(icomp==jcomp)
          {
            // diffusion
            mat(_ivar,_ivar)[count] += diff*fem.laplace(ii,jj);
          }
          count++;
        }
      }
    }
  }
}

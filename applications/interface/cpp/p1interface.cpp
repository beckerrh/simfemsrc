#include  "p1interface.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Alat/vectoronevariable.hpp"
#include  "Alat/matrixonevariable.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Solvers/meshunitwithdatainterface.hpp"
#include  <cassert>

/*--------------------------------------------------------------------------*/
P1Interface::~P1Interface() {}
P1Interface::P1Interface(): solvers::FemInterface(), _cutcells(NULL),_cutedges(NULL),_celliscut(NULL), _edgeiscut(NULL),_cutnodes(NULL),_cutcoeff(NULL),_nodesofcutcellsisin(NULL),_nodesofcutcells(NULL),_measuresofcutcells(NULL),_normalsofcutcells(NULL)
{
}
P1Interface::P1Interface( const P1Interface& P1cut): solvers::FemInterface(P1cut)
{
  assert(0);
}
P1Interface& P1Interface::operator=( const P1Interface& P1cut)
{
  assert(0);
  solvers::FemInterface::operator=(P1cut);
  return *this;
}
std::unique_ptr<solvers::FemInterface> P1Interface::clone() const
{
  return std::unique_ptr<solvers::FemInterface>(new P1Interface(*this));
}
solverEnums::fem::femtype P1Interface::getType() const {return solverEnums::fem::OWN;}
std::string P1Interface::getClassName() const
{
  return "P1Interface";
}
bool P1Interface::noIntegration() const {return true;}

/*--------------------------------------------------------------------------*/
void P1Interface::initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp)
{
  // std::cerr << "P1Interface::initFem()\n";
  solvers::FemInterface::initFem(ivar, mesh, meshinfo, ncomp);
  initCutInterface(mesh);
}


/*--------------------------------------------------------------------------*/
void P1Interface::initCutInterface(const mesh::MeshUnitInterface* mesh)
{
  // std::cerr << "P1Interface::initCutInterface()\n";
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

  nodeofcutedgelim.resize(_cutedges->size());
  nodeofcutedgeto.resize(_cutedges->size());
  nodeofcutedgelim.fill(-1);
  indicesofcutedge.resize(4,_cutedges->size());
  coefsofcutedge.resize(4,_cutedges->size());
  coefsofcutedgeweights.resize(2,_cutedges->size());

  alat::Map<int, alat::IntSet> edgesofnodes;
  alat::armaimat nodesofedges(2, _cutedges->size());
  for(int ile=0;ile<_cutedges->size();ile++)
  {
    int ie = (*_cutedges)[ile];
    int i0 = _meshinfo->nodes_of_edges(0,ie);
    int i1 = _meshinfo->nodes_of_edges(1,ie);
    edgesofnodes[i0].insert(ile);
    edgesofnodes[i1].insert(ile);
    nodesofedges(0, ile) = i0;
    nodesofedges(1, ile) = i1;
  }

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
    // if(edgesofnodes[i0].size()==1)
    // {
    //   nodeofcutedgelim[ile] = _meshinfo->nnodes+ic0;
    //   nodeofcutedgeto[ile] = _meshinfo->nnodes+ic1;
    //   coefsofcutedgeweights(0, ile) = d0;
    //   coefsofcutedgeweights(1, ile) = -d1;
    // }
    // else
    {
      nodeofcutedgelim[ile] = _meshinfo->nnodes+ic1;
      nodeofcutedgeto[ile] = _meshinfo->nnodes+ic0;
      coefsofcutedgeweights(0, ile) = d1;
      coefsofcutedgeweights(1, ile) = -d0;
    }
    coefsofcutedge(0, ile) = d0;
    coefsofcutedge(1, ile) = d1;
    coefsofcutedge(2, ile) = -d0;
    coefsofcutedge(3, ile) = -d1;
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
  // std::cerr << "nodeofcutedgelim="<<nodeofcutedgelim.t();
  // std::cerr << "indicesofcutedge="<<indicesofcutedge;
  // std::cerr << "coefsofcutedge="<<coefsofcutedge;
}

/*--------------------------------------------------------------------------*/
int P1Interface::getN() const
{
  return _cutedges->size();
}
/*--------------------------------------------------------------------------*/
int P1Interface::getNPerCell(int iK) const
{
  return 0;
}
/*--------------------------------------------------------------------------*/
void P1Interface::indicesOfCell(int iK, alat::armaivec& indices) const
{
  indices.set_size(0);
}
/*--------------------------------------------------------------------------*/
void P1Interface::setCell(int iK)
{
}

#include  "Alat/sparsitypattern.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Mesh/meshunit.hpp"
#include  "Mesh/geometryobject.hpp"
#include  <sstream>
#include  <mpi.h>

using namespace mesh;

/*---------------------------------------------------------*/
MeshUnit::~MeshUnit() {}
MeshUnit::MeshUnit() : MeshUnitInterface(), _visitor(NULL), _init_called(false)
{
  // setDebugLevel(1);
  // std::cerr << "MeshUnit::MeshUnit()\n";
}
MeshUnit::MeshUnit(const MeshUnit& meshunit) :  MeshUnitInterface(meshunit)
{
  // std::cerr << "MeshUnit::MeshUnit(const MeshUnit& meshunit)\n";
  _init_called = meshunit._init_called;
  _visitor = meshunit._visitor->clone();
  n_dim = meshunit.n_dim;
  n_nodes_per_cell = meshunit.n_nodes_per_cell;
  n_edges_per_cell = meshunit.n_edges_per_cell;
  n_sides_per_cell = meshunit.n_sides_per_cell;
  n_nodes_per_side = meshunit.n_nodes_per_side;
  _nodesandnodesofcells = meshunit._nodesandnodesofcells;
  _sidesandcells = meshunit._sidesandcells;
  _edgesandcells = meshunit._edgesandcells;
  for(BoundaryInformationMap::const_iterator p=meshunit._boundaryinformationmap.begin(); p!=meshunit._boundaryinformationmap.end(); p++)
  {
    _boundaryinformationmap[p->first] = p->second;
  }
  for(GeometryObjects::const_iterator p=meshunit._geometryobjects.begin(); p!=meshunit._geometryobjects.end(); p++)
  {
    _geometryobjects[p->first] = p->second->clone();
  }
}
void MeshUnit::init(std::shared_ptr<const mesh::MeshVisitorInterface> visitor)
{
  _visitor = visitor;
  _init_called=true;
  // std::cerr << "MeshUnit::init(std::shared_ptr<const mesh::MeshVisitorInterface> visitor)\n";
  visitor->set_dimensions(n_dim, n_nodes_per_cell, n_edges_per_cell, n_sides_per_cell, n_nodes_per_side);
  // std::cerr<< "MeshUnit:" << getInfo()<<"\n";
}
MeshUnit& MeshUnit::operator=(const MeshUnit& meshbase)
{
  assert(0);
  return *this;
}
std::string MeshUnit::getClassName() const
{
  return "mesh::MeshUnit";
}
std::string MeshUnit::getInfo() const
{
  std::stringstream ss;
  ss << "dim_nnpc_nspc_nepc_nnps: " << getDimension() << "_" << getNNodesPerCell()<< "_" << getNSidesPerCell()<<"_"<<getNEdgesPerCell()<<"_"<< getNNodesPerSide()<<" | ";
  ss << "nn_nc_ns_ne: " << getNNodes()<<"_"<< getNCells()<<"_"<< getNSides()<<"_"<< getNEdges();
  ss << "\t" << _visitor->getClassName();
  ss << "\t" << " GeometryObjects: ";
  for(GeometryObjects::const_iterator p = _geometryobjects.begin();p!=_geometryobjects.end();p++)
  {
    ss << p->second->getClassName() << " ";
  }
  return ss.str();
}

/*---------------------------------------------------------*/
alat::SparsityPattern MeshUnit::getSizes() const
{
  // 0 : NodesAndNodesOfCells
  // 1: SidesAndCells
  // 2: BoundaryInformationMap colors
  // 3: GeometryObjects enums
  int n = 5 + _boundaryinformationmap.size() + _geometryobjects.size();
  std::vector<alat::armaivec> estalavida(n);
  estalavida[0] = _nodesandnodesofcells.getSizes();
  estalavida[1] = _sidesandcells.getSizes();
  estalavida[2] = _edgesandcells.getSizes();
  alat::IntSet bdrycolors = _boundaryinformationmap.keys();
  estalavida[3] = alat::armaivec(bdrycolors.size());
  std::copy(bdrycolors.begin(), bdrycolors.end(), estalavida[2].begin());
  alat::Set<meshEnums::geomobjtype> geomenums = _geometryobjects.keys();
  estalavida[4] = alat::armaivec(geomenums.size());
  int count=0;
  for(alat::Set<meshEnums::geomobjtype>::const_iterator p=geomenums.begin();p!=geomenums.end();p++)
  {
    estalavida[4][count++] = *p;
  }
  count=0;
  for(alat::IntSet::const_iterator p=bdrycolors.begin();p!=bdrycolors.end();p++)
  {
    estalavida[5+count++] = _boundaryinformationmap[*p].getSizes();
  }
  for(alat::Set<meshEnums::geomobjtype>::const_iterator p=geomenums.begin();p!=geomenums.end();p++)
  {
    estalavida[5+count++] = _geometryobjects[*p]->getSizes();
  }
  // std::cerr << "MeshUnit::getSizes() sparsitypattern=" << estalavida[2] << "\n";
  return alat::SparsityPattern(estalavida);
}
void MeshUnit::setSizes(const alat::SparsityPattern& sparsitypattern)
{
  // std::cerr << "MeshUnit::setSizes() sparsitypattern=" << sparsitypattern << "\n";
  // _nodesandnodesofcells.setSizes(sparsitypattern.col(), sparsitypattern.rowstart(0));
  alat::armaivec::const_iterator sizes = sparsitypattern.col().begin() + sparsitypattern.rowstart(0);
  _nodesandnodesofcells.setSizes(sizes);
  sizes = sparsitypattern.col().begin() + sparsitypattern.rowstart(1);
  _sidesandcells.setSizes(sizes);
  sizes = sparsitypattern.col().begin() + sparsitypattern.rowstart(2);
  _edgesandcells.setSizes(sizes);
  alat::armaivec::const_iterator bisizesA = sparsitypattern.col().begin() + sparsitypattern.rowstart(3);
  alat::armaivec::const_iterator bisizesE = sparsitypattern.col().begin() + sparsitypattern.rowstop(3);
  int count=0;
  for(alat::armaivec::const_iterator p=bisizesA;p!=bisizesE;p++)
  {
    // std::cerr << "setting size for " << *p << "\n";
    sizes = sparsitypattern.col().begin() + sparsitypattern.rowstart(5+count++);
    _boundaryinformationmap[*p].setSizes(sizes);
  }
  alat::armaivec::const_iterator geomsizesA = sparsitypattern.col().begin() + sparsitypattern.rowstart(4);
  alat::armaivec::const_iterator geomsizesE = sparsitypattern.col().begin() + sparsitypattern.rowstop(4);
  for(alat::armaivec::const_iterator p=geomsizesA;p!=geomsizesE;p++)
  {
    meshEnums::geomobjtype type = static_cast<meshEnums::geomobjtype>(*p);
    _geometryobjects[type] = _geometryobjects_constructor.newGeometryObject(type);
    sizes = sparsitypattern.col().begin() + sparsitypattern.rowstart(5+count++);
    _geometryobjects[type]->setSizes(sizes);
  }
  // std::cerr << "BI " << sparsitypattern.rowsize(3) << "\n";
  // std::copy(sizesA, sizesE, std::ostream_iterator<int>(std::cerr, " "));
}
void MeshUnit::sendRecv(std::shared_ptr<mesh::MeshUnitInterface> meshunit, int neighbor) const
{
  int tag = 0;
  _nodesandnodesofcells.send(neighbor, tag);
  NodesAndNodesOfCells& nodesandnodesofcells_recieve = meshunit->getNodesAndNodesOfCells();
  nodesandnodesofcells_recieve.recv(neighbor, tag);

  tag = 1;
  _sidesandcells.send(neighbor, tag);
  SidesAndCells& sidesandcells_recieve = meshunit->getSidesAndCells();
  sidesandcells_recieve.recv(neighbor, tag);

  tag = 2;
  _edgesandcells.send(neighbor, tag);
  EdgesAndCells& edgesandcells_recieve = meshunit->getEdgesAndCells();
  edgesandcells_recieve.recv(neighbor, tag);

  tag = 3;
  // std::cerr << "neighbor = " << neighbor << "\n";
  // std::cerr << "_boundaryinformationmap.keys() = " << _boundaryinformationmap.keys() << "\n";
  // std::cerr << "boundaryinformationmap_recieve.keys() = " << boundaryinformationmap_recieve.keys() << "\n";
  for(BoundaryInformationMap::const_iterator p = _boundaryinformationmap.begin();p!=_boundaryinformationmap.end();p++)
  {
    p->second.send(neighbor, tag);
  }
  BoundaryInformationMap& boundaryinformationmap_recieve = meshunit->getBoundaryInformationMap();
  for(BoundaryInformationMap::iterator p = boundaryinformationmap_recieve.begin();p!=boundaryinformationmap_recieve.end();p++)
  {
    p->second.recv(neighbor, tag);
  }

  tag = 4;
  for(GeometryObjects::const_iterator p = _geometryobjects.begin();p!=_geometryobjects.end();p++)
  {
    p->second->send(neighbor, tag);
  }
  GeometryObjects& geometryobjects_recieve = meshunit->getGeometryObjects();
  for(GeometryObjects::iterator p = geometryobjects_recieve.begin();p!=geometryobjects_recieve.end();p++)
  {
    p->second->recv(neighbor, tag);
  }
}

/*---------------------------------------------------------*/
void MeshUnit::loadH5(const arma::hdf5_name& spec)
{
  // assert(_dimensions_set);
  _nodesandnodesofcells.loadH5(arma::hdf5_name(spec.filename, spec.dsname+"/"+_nodesandnodesofcells.getClassName(), spec.opts));
  _sidesandcells.loadH5(arma::hdf5_name(spec.filename, spec.dsname+"/"+_sidesandcells.getClassName(), spec.opts));
  _edgesandcells.loadH5(arma::hdf5_name(spec.filename, spec.dsname+"/"+_edgesandcells.getClassName(), spec.opts));
  //GeometryObjects
  alat::armaivec itypes;
  itypes.load(arma::hdf5_name(spec.filename,spec.dsname+"/GeometryObjects/enums", spec.opts));
  for(int i=0;i<itypes.size();i++)
  {
    meshEnums::geomobjtype type = static_cast<meshEnums::geomobjtype>(itypes[i]);
    _geometryobjects[type] = _geometryobjects_constructor.newGeometryObject(type);
    std::string name = meshEnums::geomObjTypeToString(type);
    std::string groupname = "/GeometryObjects/" + name;
    _geometryobjects[type]->loadH5(arma::hdf5_name(spec.filename, spec.dsname+groupname, spec.opts));
  }
  //BoundaryInfo
  alat::armaivec icolors;
  icolors.load(arma::hdf5_name(spec.filename,spec.dsname+"/BoundaryInformation/colors", spec.opts));
  for(int i=0;i<icolors.size();i++)
  {
    int color = icolors[i];
    std::stringstream ss;
    ss << "/BoundaryInformation/Boundary" << color;
    _boundaryinformationmap[color].loadH5(arma::hdf5_name(spec.filename, spec.dsname+ss.str(), spec.opts));
  }
}

/*---------------------------------------------------------*/
void MeshUnit::saveH5(const arma::hdf5_name& spec) const
{
  _nodesandnodesofcells.saveH5(arma::hdf5_name(spec.filename, spec.dsname+"/"+_nodesandnodesofcells.getClassName(), spec.opts));
  _sidesandcells.saveH5(arma::hdf5_name(spec.filename, spec.dsname+"/"+_sidesandcells.getClassName(), spec.opts));
  _edgesandcells.saveH5(arma::hdf5_name(spec.filename, spec.dsname+"/"+_edgesandcells.getClassName(), spec.opts));
  //GeometryObjects
  int nobj = _geometryobjects.size();
  alat::armaivec itypes(nobj);
  int count=0;
  for(GeometryObjects::const_iterator p = _geometryobjects.begin();p!=_geometryobjects.end();p++)
  {
    itypes[count++] = p->first;
  }
  itypes.save(arma::hdf5_name(spec.filename,spec.dsname+"/GeometryObjects/enums",spec.opts));
  for(GeometryObjects::const_iterator p = _geometryobjects.begin();p!=_geometryobjects.end();p++)
  {
    std::string groupname = "/GeometryObjects/" + meshEnums::geomObjTypeToString(p->first);
    p->second->saveH5(arma::hdf5_name(spec.filename, spec.dsname+groupname, spec.opts));
  }
  //BoundaryInformation
  alat::IntSet colors = _boundaryinformationmap.keys();
  alat::armaivec icolors(colors.size());
  std::copy(colors.begin(), colors.end(), icolors.begin());
  icolors.save(arma::hdf5_name(spec.filename, spec.dsname+"/BoundaryInformation/colors", spec.opts));
  for(BoundaryInformationMap::const_iterator p = _boundaryinformationmap.begin();p!=_boundaryinformationmap.end();p++)
  {
    std::stringstream ss;
    ss << "/BoundaryInformation/Boundary" << p->first;
    p->second.saveH5(arma::hdf5_name(spec.filename, spec.dsname+ss.str(), spec.opts));
  }
}


/*---------------------------------------------------------*/
const mesh::MeshVisitorInterface* MeshUnit::getVisitor() const {return _visitor.get();}

/*---------------------------------------------------------*/
int MeshUnit::getVtkCellType() const
{
  if(n_dim==1)
  {
    return 3;
  }
  else if(n_dim==2)
  {
    if(n_nodes_per_cell == 3) return 5;
    return 9;
  }
  else
  {
    if(n_nodes_per_cell == 8) return 12;
    return 10;
  }
}
int MeshUnit::getVtkSideType() const
{
  if(n_dim==1)
  {
    return 1;
  }
  else if(n_dim==2)
  {
    return 3;
  }
  else
  {
    if(n_nodes_per_cell == 8) return 9;
    return 5;
  }
}
int MeshUnit::getGmshCellType() const
{
  if(n_dim==1)
  {
    return 1;
  }
  else if(n_dim==2)
  {
    if(n_nodes_per_cell == 3) return 2;
    return 3;
  }
  else
  {
    if(n_nodes_per_cell == 8) return 5;
    return 4;
  }
}
int MeshUnit::getGmshSideType() const
{
  if(n_dim==1)
  {
    return 15;
  }
  else if(n_dim==2)
  {
    return 1;
  }
  else
  {
    if(n_nodes_per_cell == 8) return 3;
    return 2;
  }
}
int MeshUnit::getGmshSideSideType() const
{
  if(n_dim==1)
  {
    return -1;
  }
  else if(n_dim==2)
  {
    return 15;
  }
  else
  {
    return 1;
  }
}
int MeshUnit::getGmshSideSideSideType() const
{
  if(n_dim==1)
  {
    return -1;
  }
  else if(n_dim==2)
  {
    return -1;
  }
  else
  {
    return 15;
  }
}
std::string MeshUnit::getXdmfTopoType() const
{
  if(n_dim==1)
  {
    return "Polyline";
  }
  else if(n_dim==2)
  {
    if(n_nodes_per_cell==3)
    {
      return "Triangle";
    }
    else if(n_nodes_per_cell==4)
    {
      return "Quadrilateral";
    }
    else
    {
      _error_string("getXdmfTopoType","unknown n_nodes_per_cell",n_nodes_per_cell);
      return "none";
    }
  }
  else if(n_dim==3)
  {
    if(n_nodes_per_cell==4)
    {
      return "Tetrahedron";
    }
    else if(n_nodes_per_cell==8)
    {
      return "Hexahedron";
    }
    else
    {
      _error_string("getXdmfTopoType","unknown n_nodes_per_cell",n_nodes_per_cell);
      return "none";
    }
  }
  _error_string("getXdmfTopoType","unknown n_dim",n_dim);
  return "none";
}

/*---------------------------------------------------------*/
alat::armaivec MeshUnit::_getSideOfCell(int i, int ii) const
{
  // std::cerr << "n_nodes_per_cell="<<n_nodes_per_cell << " n_nodes_per_side=" << n_nodes_per_side << " n_dim="<<n_dim<<"\n";
  assert(n_nodes_per_cell==n_nodes_per_side+1);
  alat::armaivec side(n_nodes_per_side);
  for(int jj = 0; jj < n_nodes_per_cell; jj++)
  {
    side[jj] = _nodesandnodesofcells._nodes_of_cells((ii+1+jj)%n_nodes_per_cell, i);
  }
  return side;
}
alat::armaivec MeshUnit::_getEdgeOfCell(int i, int ii) const
{
  alat::armaivec side(2);
  if(n_dim==1)
  {
    for(int jj = 0; jj < 2; jj++)
    {
      side[jj] = _nodesandnodesofcells._nodes_of_cells(jj, i);
    }
  }
  else if(n_dim==2)
  {
    int iS = _sidesandcells._sides_of_cells(ii,i);
    for(int jj = 0; jj < 2; jj++)
    {
      side[jj] = _sidesandcells._nodes_of_sides(jj, iS);
    }
  }
  else if(n_dim==3)
  {
    int i1 = (2*ii)/5;
    int i2 = ii-i1+1-ii/3;
    side[0] = _nodesandnodesofcells._nodes_of_cells(i1, i);
    side[1] = _nodesandnodesofcells._nodes_of_cells(i2, i);
  }
  else
  {
    assert(0);
  }
  return side;
}

/*---------------------------------------------------------*/
const MeshUnitInterface::BoundaryInformationMap& MeshUnit::getBoundaryInformationMap()const {return _boundaryinformationmap;}
MeshUnitInterface::BoundaryInformationMap& MeshUnit::getBoundaryInformationMap() {return _boundaryinformationmap;}
const BoundaryInformation& MeshUnit::getBoundaryInformation(int color) const {return _boundaryinformationmap[color];}
BoundaryInformation& MeshUnit::getBoundaryInformation(int color) {return _boundaryinformationmap[color];}
const NodesAndNodesOfCells& MeshUnit::getNodesAndNodesOfCells() const{return _nodesandnodesofcells;}
NodesAndNodesOfCells& MeshUnit::getNodesAndNodesOfCells(){return _nodesandnodesofcells;}
const SidesAndCells& MeshUnit::getSidesAndCells() const {return _sidesandcells;}
SidesAndCells& MeshUnit::getSidesAndCells(){return _sidesandcells;}
const EdgesAndCells& MeshUnit::getEdgesAndCells() const{return _edgesandcells;}
EdgesAndCells& MeshUnit::getEdgesAndCells(){return _edgesandcells;}

const MeshUnitInterface::GeometryObjects& MeshUnit::getGeometryObjects() const {return _geometryobjects;}
MeshUnitInterface::GeometryObjects& MeshUnit::getGeometryObjects() {return _geometryobjects;}

/*---------------------------------------------------------*/
void MeshUnit::addGeometryObject(meshEnums::geomobjtype type, const GeometryConstructorInterface* geometryconstructor)
{
  // std::cerr << "MeshUnit::addGeometryObject() "<< meshEnums::geomObjTypeToString(type) << "\n";
  if(_debug_level)  std::cerr << "MeshUnit::addGeometryObject() _debug_level = "<< _debug_level << "\n";
  if(_geometryobjects.find(type)!=_geometryobjects.end())
  {
    _error_string("addGeometryObject", _geometryobjects[type]->getClassName(), "exists already");
  }
  _geometryobjects[type] = _geometryobjects_constructor.newGeometryObject(type);
  // std::cerr << "MeshUnit::addGeometryObject() construct=\n";
  if(_debug_level) _geometryobjects[type]->setDebugLevel(_debug_level);
  // std::cerr << "MeshUnit::addGeometryObject() construct=\n";
  _geometryobjects[type]->construct(this, geometryconstructor);
  // std::cerr << "MeshUnit::addGeometryObject() construct=\n";
  if(_debug_level) _geometryobjects[type]->setDebugLevel(0);
  // std::cerr << "MeshUnit::addGeometryObject() "<< meshEnums::geomObjTypeToString(type) << "\n";
}
bool MeshUnit::geometryObjectExists(meshEnums::geomobjtype type) const
{
  return _geometryobjects.hasKey(type);
}
std::shared_ptr<const GeometryObject> MeshUnit::getGeometryObject(meshEnums::geomobjtype type) const
{
  return _geometryobjects[type];
}
std::shared_ptr<GeometryObject> MeshUnit::getGeometryObject(meshEnums::geomobjtype type)
{
  return _geometryobjects[type];
}

/*---------------------------------------------------------*/
int MeshUnit::getDimension() const
{
  return n_dim;
}
int MeshUnit::getNNodesPerCell() const
{
  return n_nodes_per_cell;
}
/*---------------------------------------------------------*/
int MeshUnit::getNNodesPerSide() const
{
  return n_nodes_per_side;
}
/*---------------------------------------------------------*/
int MeshUnit::getNSidesPerCell() const
{
  return n_sides_per_cell;
}
/*---------------------------------------------------------*/
int MeshUnit::getNEdgesPerCell() const
{
  return n_edges_per_cell;
}

/*---------------------------------------------------------*/
int MeshUnit::getNNodes() const
{
  return _nodesandnodesofcells._nodes.n_cols;
}
int MeshUnit::getNCells() const
{
  return  _nodesandnodesofcells._nodes_of_cells.n_cols;
}
int MeshUnit::getNSides() const
{
  return MeshUnit::_sidesandcells._nodes_of_sides.n_cols;
}
int MeshUnit::getNEdges() const
{
  if(n_dim==1) {return getNCells();}
  else if(n_dim==2) {return getNSides();}
  return MeshUnit::_edgesandcells._nodes_of_edges.n_cols;
}

/*---------------------------------------------------------*/
const arma::mat& MeshUnit::getNodes() const
{
  return _nodesandnodesofcells._nodes;
}
arma::mat& MeshUnit::getNodes()
{
  return _nodesandnodesofcells._nodes;
}
/*---------------------------------------------------------*/
alat::armaimat& MeshUnit::getCells()
{
  return _nodesandnodesofcells._nodes_of_cells;
}
const alat::armaimat& MeshUnit::getCells() const
{
  return _nodesandnodesofcells._nodes_of_cells;
}
/*---------------------------------------------------------*/
alat::armaimat& MeshUnit::getSides()
{
  return _sidesandcells._nodes_of_sides;
}
const alat::armaimat& MeshUnit::getSides() const
{
  return _sidesandcells._nodes_of_sides;
}
/*---------------------------------------------------------*/
alat::armaimat& MeshUnit::getSidesOfCells()
{
  return _sidesandcells._sides_of_cells;
}
const alat::armaimat& MeshUnit::getSidesOfCells() const
{
  return _sidesandcells._sides_of_cells;
}
/*---------------------------------------------------------*/
alat::armaimat&  MeshUnit::getCellsOfSides()
{
  return _sidesandcells._cells_of_sides;
}
const alat::armaimat& MeshUnit::getCellsOfSides() const
{
  return _sidesandcells._cells_of_sides;
}
/*---------------------------------------------------------*/
const MeshUnitInterface::Cell MeshUnit::getCell(int i) const
{
  return _nodesandnodesofcells._nodes_of_cells.col(i);
}
/*---------------------------------------------------------*/
const MeshUnitInterface::Side MeshUnit::getSide(int i) const
{
  return _sidesandcells._nodes_of_sides.col(i);
}
/*---------------------------------------------------------*/
alat::Node MeshUnit::getNodeOfCell(int iK) const
{
  alat::Node v;
  double d = 1.0/double(n_nodes_per_cell);
  for(int ii = 0; ii < n_nodes_per_cell; ii++)
  {
    int i = _nodesandnodesofcells._nodes_of_cells(ii,iK);
    v.x() += d*_nodesandnodesofcells._nodes(0,i);
    v.y() += d*_nodesandnodesofcells._nodes(1,i);
    v.z() += d*_nodesandnodesofcells._nodes(2,i);
  }
  return v;
}
alat::Node MeshUnit::getNodeOfSide(int iS) const
{
  const arma::mat& _nodes = _nodesandnodesofcells._nodes;
  alat::Node v;
  double d = 1.0/double(n_nodes_per_side);
  for(int ii = 0; ii < n_nodes_per_side; ii++)
  {
    int i = _sidesandcells._nodes_of_sides(ii,iS);
    v.x() += d*_nodesandnodesofcells._nodes(0,i);
    v.y() += d*_nodesandnodesofcells._nodes(1,i);
    v.z() += d*_nodesandnodesofcells._nodes(2,i);
  }
  return v;
}

/*--------------------------------------------------------------------------*/
void MeshUnit::writeVtk(std::string filename) const
{
  const arma::mat& _nodes = _nodesandnodesofcells._nodes;
  const alat::armaimat& _nodes_of_cells = _nodesandnodesofcells._nodes_of_cells;

  std::string name = filename;
  // name += ".vtk";

  std::ofstream file( name.c_str() );
  assert( file.is_open() );


  file<<"# vtk DataFile Version 4.0 "<<std::endl;
  file<<"output from QuadrilateralMesh"<<std::endl;
  file<<"ASCII"<<std::endl;
  //     file << "binary" << std::endl;
  file<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
  file<<std::endl;

  int nn = getNNodes();

  file<<"POINTS "<<nn;
  file<<" FLOAT"<<std::endl;
  for(int i = 0; i < nn; i++)
  {
    file<<_nodes(0,i)<<" "<<_nodes(1,i)<<" "<<_nodes(2,i)<<" "<<std::endl;
  }
  file<<std::endl;

  int ne = getNCells();
  int nle = n_nodes_per_cell;
  int lenght = ne*( nle+1 );

  file<<std::endl<<"CELLS "<<ne<<" "<<lenght<<std::endl;

  for(int ie = 0; ie < ne; ie++)
  {
    file<<nle<<" ";
    for(int ii = 0; ii < nle; ii++)
    {
      file<<_nodes_of_cells(ii, ie)<<" ";
    }
    file<<std::endl;
  }
  file<<std::endl<<"CELL_TYPES "<<ne<<std::endl;
  for(int c = 0; c < ne; c++)
  {
    file<<getVtkCellType()<<" ";
  }
  file<<std::endl;

  file.close();
}

/*--------------------------------------------------------------------------*/
void MeshUnit::writeBoundaryVtk(std::string filename) const
{
  const arma::mat& _nodes = _nodesandnodesofcells._nodes;
  const alat::armaimat& _nodes_of_cells = _nodesandnodesofcells._nodes_of_cells;
  const alat::armaimat& _nodes_of_sides = _sidesandcells._nodes_of_sides;

  std::string name = filename;
  // name += "-boundary.vtk";

  std::ofstream file( name.c_str() );
  assert( file.is_open() );

  file<<"# vtk DataFile Version 4.0 "<<std::endl;
  file<<"output from QuadrilateralMesh"<<std::endl;
  file<<"ASCII"<<std::endl;
  //     file << "binary" << std::endl;
  file<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
  file<<std::endl;

  int nn = getNNodes();

  file<<"POINTS "<<nn;
  file<<" FLOAT"<<std::endl;
  for(int i = 0; i < nn; i++)
  {
    file<<_nodes(0,i)<<" "<<_nodes(1,i)<<" "<<_nodes(2,i)<<" "<<std::endl;
  }
  file<<std::endl;

  int nsides=0;
  for(BoundaryInformationMap::const_iterator p=_boundaryinformationmap.begin();p!=_boundaryinformationmap.end();p++)
  {
    // alat::armaimat& cells_on_bdry = p->second.get()->getBoundaryMeshInformation().getCellsOnBdryOfPlain();
    nsides += p->second.size();
  }
  int lenght = nsides*( n_nodes_per_side+1 );
  file<<std::endl<<"CELLS "<<nsides<<" "<<lenght<<std::endl;

  for(BoundaryInformationMap::const_iterator p=_boundaryinformationmap.begin();p!=_boundaryinformationmap.end();p++)
  {
    const alat::armaimat& cells_on_bdry = p->second.getCellsOnBdryOfPlain();
    for(int j = 0; j < cells_on_bdry.n_cols; j++)
    {
      file<<n_nodes_per_side<<" ";
      for(int ii = 0; ii < n_nodes_per_side; ii++)
      {
        file<<_nodes_of_sides(ii,cells_on_bdry(1,j))<<" ";
      }
      file<<std::endl;
    }
  }
  file<<std::endl<<"CELL_TYPES "<<nsides<<std::endl;
  for(int c = 0; c < nsides; c++)
  {
    file<<getVtkSideType()<<" ";
  }
  file<<std::endl;
  file<<std::endl<<"CELL_DATA "<<nsides<<std::endl;
  file<<std::endl<<"SCALARS "<<" bdry_colors "<<" int "<<1<<std::endl;
  file<<std::endl<<"LOOKUP_TABLE default"<<std::endl;
  for(BoundaryInformationMap::const_iterator p=_boundaryinformationmap.begin();p!=_boundaryinformationmap.end();p++)
  {
    const alat::armaimat& cells_on_bdry = p->second.getCellsOnBdryOfPlain();
    for(int j = 0; j < cells_on_bdry.n_cols; j++)
    {
      file<<p->first<<" ";
    }
    file<<std::endl;
  }
  file.close();
}

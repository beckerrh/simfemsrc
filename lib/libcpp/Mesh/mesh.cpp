#include  "Mesh/mesh.hpp"
#include  "Mesh/sidesconstructor.hpp"
#include  "Mesh/edgesconstructor.hpp"
#include  "Mesh/readergmshnew.hpp"
#include  "Mesh/meshvisitorpoint.hpp"
#include  "Mesh/meshvisitorline.hpp"
#include  "Mesh/meshvisitortetrahedral.hpp"
#include  "Mesh/meshvisitortriangle.hpp"
#include  "Mesh/sidesconstructor.hpp"
  #include  <sstream>
#include  <mpi.h>

using namespace mesh;

/*---------------------------------------------------------*/
Mesh::~Mesh() {}
Mesh::Mesh() : mesh::MeshInterface()
{
  // std::cerr << "Mesh::Mesh()\n";
}
Mesh::Mesh(const Mesh& mesh) : mesh::MeshInterface(mesh), _plainmesh(mesh._plainmesh)
{
  // std::cerr << "Mesh::Mesh(const Mesh& meshbase)\n";
  _partion_id = mesh._partion_id;
  _construct_bdrymeshes = mesh._construct_bdrymeshes;
  for(BoundaryMeshUnitsMap::const_iterator p=mesh._boundarymeshunits.begin(); p!=mesh._boundarymeshunits.end(); p++)
  {
    _boundarymeshunits[p->first] = p->second;
  }
  for(InterfaceMeshUnitsMap::const_iterator p=mesh._interfacemeshunits.begin(); p!=mesh._interfacemeshunits.end(); p++)
  {
    _interfacemeshunits[p->first] = p->second;
  }
}
/*---------------------------------------------------------*/
void Mesh::exchangeInterfaces(int nproc)
{
  alat::IntSet keys = _interfacemeshunits.keys();
  for(alat::IntSet::const_iterator p=keys.begin(); p!=keys.end(); p++)
  {
    int neighbor = -*p;
    _interfacemeshunits[neighbor] = std::unique_ptr<InterfaceMeshUnit>(new InterfaceMeshUnit());
    _interfacemeshunits[neighbor]->init(_plainmesh.getVisitor()->clone(), &_plainmesh);
  }
  _interfacemeshunits.exchangeInterfaces(_partion_id, nproc);
}

/*---------------------------------------------------------*/
std::unique_ptr<MeshInterface> Mesh::create(std::string name, int partion_id, bool construct_bdrymeshes)
{
  return create(meshEnums::stringTomeshType(name), partion_id, construct_bdrymeshes);
}
std::unique_ptr<MeshInterface> Mesh::create(meshEnums::meshtype type, int partion_id, bool construct_bdrymeshes)
{
  Mesh* newmesh = new Mesh;
  if(type==meshEnums::PointMesh)
  {
    newmesh->init(std::unique_ptr<MeshVisitorInterface>(new MeshVisitorPoint), partion_id, construct_bdrymeshes);
  }
  else if(type==meshEnums::LineMesh)
  {
    newmesh->init(std::unique_ptr<MeshVisitorInterface>(new MeshVisitorLine), partion_id, construct_bdrymeshes);
  }
  else if(type==meshEnums::TriangleMesh)
  {
    newmesh->init(std::unique_ptr<MeshVisitorInterface>(new MeshVisitorTriangle), partion_id, construct_bdrymeshes);
  }
  else if(type==meshEnums::TetrahedralMesh)
  {
    newmesh->init(std::unique_ptr<MeshVisitorInterface>(new MeshVisitorTetrahedral), partion_id, construct_bdrymeshes);
  }
  else
  {
    newmesh->_error_string("create", "unknown mesh type=", meshEnums::meshTypeToString(type));
  }
  return std::unique_ptr<MeshInterface>(newmesh);
}

Mesh& Mesh::operator=(const Mesh& meshbase)
{
  assert(0);
  return *this;
}
void Mesh::init(std::shared_ptr<const mesh::MeshVisitorInterface> visitor, int partion_id, bool construct_bdrymeshes)
{
  _partion_id = partion_id;
  _construct_bdrymeshes = construct_bdrymeshes;
  _plainmesh.init(visitor);
}
std::string Mesh::getClassName() const
{
  return "mesh::Mesh";
}
/*---------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
std::ostream& mesh::operator<<(std::ostream &os, const Mesh& mesh)
{
	os << mesh.getInfo();
	return os;
}
int Mesh::getPartionId() const {return _partion_id;}
void Mesh::setPartionId(int id) {_partion_id=id;}
void Mesh::addGeometryObjectByName(std::string name)
{
	addGeometryObject(meshEnums::stringTogeomObjType(name));
}
void Mesh::addGeometryObject(meshEnums::geomobjtype type, const GeometryConstructorInterface* geometryconstructor)
{
  _plainmesh.addGeometryObject(type, geometryconstructor);
  if(_construct_bdrymeshes==false) return;
  if(_plainmesh.getDimension()==1) return;
  for(BoundaryMeshUnitsMap::const_iterator p=_boundarymeshunits.begin(); p!=_boundarymeshunits.end();p++)
  {
    // std::cerr << "addGeometryObject() "<< meshEnums::geomObjTypeToString(type)<< " for color " << p->first << "\n";
    p->second->addGeometryObject(type, geometryconstructor);
  }
  for(InterfaceMeshUnitsMap::const_iterator p=_interfacemeshunits.begin(); p!=_interfacemeshunits.end();p++)
  {
    // std::cerr << "addGeometryObject() "<< meshEnums::geomObjTypeToString(type)<< " for color " << p->first << "\n";
    // p->second->setDebugLevel(2);
    p->second->addGeometryObject(type, geometryconstructor);
    // p->second->setDebugLevel(0);
  }
}
/*---------------------------------------------------------*/
int Mesh::getNCells() const
{
  return _plainmesh.getNCells();
}
/*---------------------------------------------------------*/
int Mesh::getDimension() const
{
  return _plainmesh.getDimension();
}
/*---------------------------------------------------------*/
std::string Mesh::getInfo() const
{
  std::string sep("-----------------\n");
  std::string info = sep + "PlainMesh: " + _plainmesh.getInfo() + "\n";
  if(_construct_bdrymeshes)
  {
    for(BoundaryMeshUnitsMap::const_iterator p=_boundarymeshunits.begin(); p!=_boundarymeshunits.end();p++)
    {
      std::stringstream ss;
      ss << "BoundaryMesh " << p->first << " : " << p->second->getInfo() << "\n";
      info += ss.str();
    }
    for(InterfaceMeshUnitsMap::const_iterator p=_interfacemeshunits.begin(); p!=_interfacemeshunits.end();p++)
    {
      std::stringstream ss;
      ss << "InterfaceMesh " << p->first << " : " << p->second->getInfo() << "\n";
      info += ss.str();
    }
  }
  return info+ sep;
}

/*---------------------------------------------------------*/
const MeshUnitInterface* Mesh::getPlainMesh() const {return &_plainmesh;}
const MeshUnitInterface* Mesh::getBoundaryMesh(int color) const {return _boundarymeshunits[color].get();}
const MeshUnitInterface* Mesh::getInterfaceMesh(int color) const {return _interfacemeshunits[color].get();}
const mesh::MeshInterface::BoundaryMeshUnitsMap& Mesh::getBoundaryMeshUnitsMap() const {return _boundarymeshunits;}
const InterfaceMeshUnitsMap& Mesh::getInterfaceMeshUnitsMap() const {return _interfacemeshunits;}
MeshUnitInterface* Mesh::getPlainMesh(){return &_plainmesh;}
mesh::MeshInterface::BoundaryMeshUnitsMap& Mesh::getBoundaryMeshUnitsMap(){return _boundarymeshunits;}
InterfaceMeshUnitsMap& Mesh::getInterfaceMeshUnitsMap(){return _interfacemeshunits;}

/*---------------------------------------------------------*/
void Mesh::loadH5(std::string filename)
{
  arma::hdf5_name spec(filename,"PlainMesh");
  _plainmesh.loadH5(spec);
  //BoundaryMeshes
  alat::armaivec icolors;
  icolors.load(arma::hdf5_name(filename,"BoundaryMeshes/colors", spec.opts));
  for(int i=0;i<icolors.size();i++)
  {
    int color = icolors[i];
    std::stringstream ss;
    ss << "BoundaryMeshes/Mesh" << color;
    // _boundarymeshunits[color] = newBoundaryMesh(this);
    _boundarymeshunits[color] = std::unique_ptr<BoundaryMeshUnit>(new BoundaryMeshUnit());
    _boundarymeshunits[color]->init(_plainmesh.getVisitor()->newBoundaryVisitor(), &_plainmesh);
    _boundarymeshunits[color]->loadH5(arma::hdf5_name(filename, ss.str(), spec.opts));
  }
  //InterfaceMeshes
  icolors.load(arma::hdf5_name(filename,"InterfaceMeshes/colors", spec.opts));
  for(int i=0;i<icolors.size();i++)
  {
    int color = icolors[i];
    std::stringstream ss;
    ss << "InterfaceMeshes/Mesh" << color;
    _interfacemeshunits[color] = std::unique_ptr<InterfaceMeshUnit>(new InterfaceMeshUnit());
    _interfacemeshunits[color]->init(_plainmesh.getVisitor()->clone(), &_plainmesh);
    _interfacemeshunits[color]->loadH5(arma::hdf5_name(filename, ss.str(), spec.opts));
  }
}

/*---------------------------------------------------------*/
void Mesh::saveH5(std::string filename) const
{
  #include <stdio.h>
  remove( filename.c_str() );
  arma::hdf5_name spec(filename,"PlainMesh",arma::hdf5_opts::append);
  _plainmesh.saveH5(spec);
  //BoundaryMeshes
  alat::IntSet colors = _boundarymeshunits.keys();
  alat::armaivec icolors(colors.size());
  std::copy(colors.begin(), colors.end(), icolors.begin());
  // std::cerr << "icolors=" << icolors.t();
  if(not _construct_bdrymeshes)
  {
    icolors.set_size(0);
  }
  // std::cerr << "_boundarymeshunits.size()=" << _boundarymeshunits.size()<<"\n";
  assert(icolors.size()==_boundarymeshunits.size());
  icolors.save(arma::hdf5_name(spec.filename, "BoundaryMeshes/colors", spec.opts));
  for(BoundaryMeshUnitsMap::const_iterator p = _boundarymeshunits.begin();p!=_boundarymeshunits.end();p++)
  {
    std::stringstream ss;
    ss << "BoundaryMeshes/Mesh" << p->first;
    // std::cerr << "### saveH5() bdry mesh " << spec.dsname+ss.str() << "\n";
    p->second->saveH5(arma::hdf5_name(spec.filename, ss.str(), spec.opts));
  }
  //InterfaceMeshes
  colors = _interfacemeshunits.keys();
  icolors.set_size(colors.size());
  std::copy(colors.begin(), colors.end(), icolors.begin());
  icolors.save(arma::hdf5_name(spec.filename, "InterfaceMeshes/colors", spec.opts));
  for(InterfaceMeshUnitsMap::const_iterator p = _interfacemeshunits.begin();p!=_interfacemeshunits.end();p++)
  {
    std::stringstream ss;
    ss << "InterfaceMeshes/Mesh" << p->first;
    p->second->saveH5(arma::hdf5_name(spec.filename, ss.str(), spec.opts));
  }
}

/*--------------------------------------------------------------------------*/
void Mesh::writeVtk(std::string filename) const
{
  _plainmesh.writeVtk(filename);
}

/*--------------------------------------------------------------------------*/
void Mesh::writeBoundaryVtk(std::string filename) const
{
  _plainmesh.writeBoundaryVtk(filename);
}

/*--------------------------------------------------------------------------*/
void Mesh::_callSidesConstructor(MeshUnitInterface* mesh, const MeshUnitInterface::BoundarySideToColorForSidesConstructor& bsides, const alat::Map<int, alat::armaivec>* partition_to_cells)
{
  int nnodes_per_side = mesh->getNNodesPerSide();
  if(nnodes_per_side==1)
  {
    SidesConstructor<1> sidesconstructor(mesh);
    sidesconstructor.constructSidesFromCells(bsides, partition_to_cells);
  }
  else if(nnodes_per_side==2)
  {
    SidesConstructor<2> sidesconstructor(mesh);
    sidesconstructor.constructSidesFromCells(bsides, partition_to_cells);
  }
  else if(nnodes_per_side==3)
  {
    SidesConstructor<3> sidesconstructor(mesh);
    sidesconstructor.constructSidesFromCells(bsides, partition_to_cells);
  }
  else
  {
    _error_string("_callSidesConstructor", "unknown nnodes_per_side", nnodes_per_side);
  }
  EdgesConstructor edgesconstructor(mesh);
  bool debug=false;
  edgesconstructor.constructEdgesFromCells(debug);
}

/*---------------------------------------------------------*/
void Mesh::readGmsh(std::string filename)
{
  ReaderGmshNew readergmsh;
  readergmsh.read(filename);
  // std::cerr << "mesh id =" << this->getPartionId() << "\n";
  readergmsh.setMesh(&_plainmesh);
  // std::cerr << "After ReaderGmsh\n";
  // std::cerr << "readergmsh.color_to_bdry1 " << readergmsh.color_to_bdry1 << "\n";
  // std::cerr << "readergmsh.color_to_bdry2 " << readergmsh.color_to_bdry2 << "\n";
  // std::cerr << "readergmsh.color_to_bdry3 " << readergmsh.color_to_bdry3 << "\n";
  // std::cerr << "readergmsh.partition_to_cells " << readergmsh.partition_to_cells << "\n";
  _callSidesConstructor(&_plainmesh, readergmsh.color_to_bdry1, &readergmsh.partition_to_cells);
  const MeshUnitInterface::BoundaryInformationMap& bimap = _plainmesh.getBoundaryInformationMap();
  // std::cerr << "_plainmesh.getBoundaryInformationMap() " << bimap << "\n";
  for(MeshUnitInterface::BoundaryInformationMap::const_iterator p=bimap.begin();p!=bimap.end();p++)
  {
    if(p->first<0)
    {
      // std::cerr << "readGmsh() "  << this->getPartionId() << " Constructing interface mesh for color " << p->first << "\n";
      // std::cerr << "readGmsh() BoundaryInformation " << p->second << "\n";
      _interfacemeshunits[p->first] = std::unique_ptr<InterfaceMeshUnit>(new InterfaceMeshUnit());
      _interfacemeshunits[p->first]->init(_plainmesh.getVisitor()->clone(), &_plainmesh);
      constructInterfacemesh(*_interfacemeshunits[p->first], p->second, readergmsh.color_to_bdry1);
    }
    // else
    {
      if(_construct_bdrymeshes)
      {
        // std::cerr << "readGmsh() Constructing boundary mesh for color " << p->first << "\n";
        // std::cerr << "readGmsh() BoundaryInformation " << p->second << "\n";
        _boundarymeshunits[p->first] = std::unique_ptr<BoundaryMeshUnit>(new BoundaryMeshUnit());
        _boundarymeshunits[p->first]->init(_plainmesh.getVisitor()->newBoundaryVisitor(), &_plainmesh);
        constructBdrymesh(*_boundarymeshunits[p->first], p->second, readergmsh.color_to_bdry2);
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
void Mesh::constructInterfacemesh(InterfaceMeshUnit& mesh, const BoundaryInformation& boundaryinformation, const alat::Map<int, alat::armaimat>& color_to_bdry1)
{
  // std::cerr <<"@Me " << mesh.getInfo()<<"\n";
  // std::cerr <<"@parent " << parentmesh->getInfo()<<"\n";
  const MeshUnitInterface* plainmesh = mesh.getPlainMesh();
  InterfaceMeshInfo& interfacemeshinfo = mesh.getInterfaceMeshInfo();

  // etape 1 : construct nodes
  arma::mat& nodes = mesh.getNodesAndNodesOfCells().getNodes();
  alat::armaivec&  imi_nodes_to_plain = interfacemeshinfo.getNodesToPlain();
  const alat::armaimat& cells_on_bdry = boundaryinformation.getCellsOnBdryOfPlain();
  const alat::armaimat& nodes_of_cells_plain = plainmesh->getNodesAndNodesOfCells().getNodesOfCells();
  const alat::armaimat& sides_of_cells_plain = plainmesh->getSidesAndCells().getSidesOfCells();
  const arma::mat& nodes_plain = plainmesh->getNodesAndNodesOfCells().getNodes();
  int n_nodes_per_cell_plain = plainmesh->getNNodesPerCell();
  int n_sides_per_cell_plain = plainmesh->getNSidesPerCell();
  assert(n_nodes_per_cell_plain==mesh.getNNodesPerCell());

  alat::IntMap nodes_global_to_local, sides_global_to_local;
  int count_nodes=0, count_sides=0;
  for(int i = 0; i < cells_on_bdry.n_cols; i++)
  {
    int iK = cells_on_bdry(0,i);
    // int iS = cells_on_bdry(1,i);
    // int iil = cells_on_bdry(2,i);
    for(int ii=0;ii<n_nodes_per_cell_plain;ii++)
    {
      int iN = nodes_of_cells_plain(ii,iK);
      alat::IntMap::const_iterator pN = nodes_global_to_local.find(iN);
      if(pN==nodes_global_to_local.end())
      {
        nodes_global_to_local[iN] = count_nodes++;
      }
    }
    for(int ii=0;ii<n_sides_per_cell_plain;ii++)
    {
      int iS = sides_of_cells_plain(ii,iK);
      alat::IntMap::const_iterator pS = sides_global_to_local.find(iS);
      if(pS==sides_global_to_local.end())
      {
        sides_global_to_local[iS] = count_sides++;
      }
    }
  }
  // std::cerr << "_nodes_global_to_local="<<_nodes_global_to_local<<"\n";
  imi_nodes_to_plain.set_size(nodes_global_to_local.size());
  nodes.set_size(3,count_nodes);
  for(alat::IntMap::const_iterator p= nodes_global_to_local.begin(); p!= nodes_global_to_local.end();p++)
  {
    nodes.col(p->second) = nodes_plain.col(p->first);
    imi_nodes_to_plain[p->second] = p->first;
  }

  // etape 2 : construct cells
  alat::armaimat& nodes_of_cells = mesh.getNodesAndNodesOfCells().getNodesOfCells();
  alat::armaivec& imi_cells_to_plain = interfacemeshinfo.getCellsToPlain();
  alat::IntMap cells_global_to_local;

  nodes_of_cells.set_size(mesh.getNNodesPerCell(), cells_on_bdry.n_cols);
  imi_cells_to_plain.set_size(cells_on_bdry.n_cols);
  for(int i = 0; i < cells_on_bdry.n_cols; i++)
  {
    int iK = cells_on_bdry(0,i);
    for(int ii=0;ii<n_nodes_per_cell_plain;ii++)
    {
      int iN = nodes_of_cells_plain(ii,iK);
      nodes_of_cells(ii, i) = nodes_global_to_local[iN];
    }
    imi_cells_to_plain[i] = iK;
    cells_global_to_local[iK] = i;
  }

  // etape 3 : construct sides
  alat::armaimat& sides_of_cells = mesh.getSidesAndCells().getSidesOfCells();
  alat::armaimat& cells_of_sides = mesh.getSidesAndCells().getCellsOfSides();
  alat::armaimat& nodes_of_sides = mesh.getSidesAndCells().getNodesOfSides();
  alat::armaivec& imi_sides_to_plain = interfacemeshinfo.getSidesToPlain();

  const alat::armaimat& cells_of_sides_plain = plainmesh->getSidesAndCells().getCellsOfSides();
  const alat::armaimat& nodes_of_sides_plain = plainmesh->getSidesAndCells().getNodesOfSides();
  int n_nodes_per_side_plain = plainmesh->getNNodesPerSide();
  assert(n_nodes_per_side_plain==mesh.getNNodesPerSide());

  imi_sides_to_plain.set_size(count_sides);
  nodes_of_sides.set_size(n_nodes_per_side_plain, count_sides);
  cells_of_sides.set_size(2, count_sides);
  sides_of_cells.set_size(n_sides_per_cell_plain, cells_on_bdry.n_cols);
  for(alat::IntMap::const_iterator p= sides_global_to_local.begin(); p!= sides_global_to_local.end();p++)
  {
    imi_sides_to_plain[p->second] = p->first;
    for(int ii=0;ii<n_nodes_per_side_plain;ii++)
    {
      int iN = nodes_of_sides_plain(ii,p->first);
      assert(nodes_global_to_local.find(iN)!=nodes_global_to_local.end());
      nodes_of_sides(ii, p->second) = nodes_global_to_local[iN];
    }
    for(int ii=0;ii<2;ii++)
    {
      int iK = cells_of_sides_plain(ii,p->first);
      if(cells_global_to_local.find(iK)!=cells_global_to_local.end())
      {
        cells_of_sides(ii,p->second) = cells_global_to_local[iK];
      }
      else
      {
        cells_of_sides(ii,p->second) = -1;
      }
    }
  }
  // attention il faut rectifier si des cells n'existe plus
  for(int i = 0; i < count_sides; i++)
  {
    if(cells_of_sides(0,i)==-1)
    {
      assert(cells_of_sides(1,i)!=-1);
      cells_of_sides(0,i) = cells_of_sides(1,i);
      cells_of_sides(1,i) = -1;
    }
  }
  for(int i = 0; i < cells_on_bdry.n_cols; i++)
  {
    int iK = cells_on_bdry(0,i);
    for(int ii=0;ii<n_sides_per_cell_plain;ii++)
    {
      int iS = sides_of_cells_plain(ii,iK);
      assert(sides_global_to_local.find(iS)!=sides_global_to_local.end());
      sides_of_cells(ii, i) = sides_global_to_local[iS];
    }
  }

  // etape 4 : construct boundaryinfo
  MeshUnitInterface::BoundaryInformationMap& bimap = mesh.getBoundaryInformationMap();
  const MeshUnitInterface::BoundaryInformationMap& bimap_plain = plainmesh->getBoundaryInformationMap();

  for(MeshUnitInterface::BoundaryInformationMap::const_iterator p=bimap_plain.begin();p!=bimap_plain.end();p++)
  {
    const alat::armaimat& _cells_on_bdry_of_plain = p->second.getCellsOnBdryOfPlain();
    int count = 0;
    for(int i = 0; i < _cells_on_bdry_of_plain.n_cols; i++)
    {
      int iK = _cells_on_bdry_of_plain(0,i);
      if(cells_global_to_local.find(iK)!=cells_global_to_local.end()) {count++;}
    }
    if(count==0) {continue;}
    alat::armaimat& _cells_on_bdry = bimap[p->first].getCellsOnBdryOfPlain();
    _cells_on_bdry.set_size(3,count);
    count=0;
    for(int i = 0; i < _cells_on_bdry_of_plain.n_cols; i++)
    {
      int iK = _cells_on_bdry_of_plain(0,i);
      int iS = _cells_on_bdry_of_plain(1,i);
      int iil = _cells_on_bdry_of_plain(2,i);
      if(cells_global_to_local.find(iK)!=cells_global_to_local.end())
      {
        _cells_on_bdry(0,count) = cells_global_to_local[iK];
        _cells_on_bdry(1,count) = sides_global_to_local[iS];
        _cells_on_bdry(2,count) = iil;
        count++;
      }
    }
  }


}

/*--------------------------------------------------------------------------*/
void Mesh::constructBdrymesh(BoundaryMeshUnit& mesh, const BoundaryInformation& boundaryinformation, const alat::Map<int, alat::armaimat>& color_to_bdry2)
{
  // std::cerr <<"@Me " << mesh.getInfo()<<"\n";
  // std::cerr <<"@parent " << parentmesh->getInfo()<<"\n";
  const MeshUnitInterface* parentmesh = mesh.getParentMesh();
  BoundaryMeshInfo& boundarymeshinfo = mesh.getBoundaryMeshInfo();
  alat::armaivec&  bmi_nodes_to_parent = boundarymeshinfo.getNodesToParent();
  alat::armaimat&  bmi_cells_to_parent = boundarymeshinfo.getCellsToParent();
  arma::mat& nodes = mesh.getNodesAndNodesOfCells().getNodes();
  alat::armaimat& nodes_of_cells = mesh.getNodesAndNodesOfCells().getNodesOfCells();
  const arma::mat& nodes_parent = parentmesh->getNodesAndNodesOfCells().getNodes();
  int n_nodes_per_side_parent = parentmesh->getNNodesPerSide();
  const alat::armaimat& nodes_of_sides_parents = parentmesh->getSidesAndCells().getNodesOfSides();
  const alat::armaimat& cells_on_bdry = boundaryinformation.getCellsOnBdryOfPlain();
  assert(n_nodes_per_side_parent==mesh.getNNodesPerCell());

  alat::IntMap nodes_global_to_local;
  int count_nodes=0;
  for(int i = 0; i < cells_on_bdry.n_cols; i++)
  {
    int iK = cells_on_bdry(0,i);
    int iS = cells_on_bdry(1,i);
    int iil = cells_on_bdry(2,i);
    for(int ii=0;ii<n_nodes_per_side_parent;ii++)
    {
      int iN = nodes_of_sides_parents(ii,iS);
      alat::IntMap::const_iterator pN = nodes_global_to_local.find(iN);
      if(pN==nodes_global_to_local.end())
      {
        nodes_global_to_local[iN] = count_nodes++;
      }
    }
  }
  // std::cerr << "_nodes_global_to_local="<<_nodes_global_to_local<<"\n";
  bmi_nodes_to_parent.set_size(nodes_global_to_local.size());
  nodes.set_size(3,count_nodes);
  for(alat::IntMap::const_iterator p= nodes_global_to_local.begin(); p!= nodes_global_to_local.end();p++)
  {
    nodes.col(p->second) = nodes_parent.col(p->first);
    bmi_nodes_to_parent[p->second] = p->first;
  }
  // std::cerr << "nodes_parent="<<nodes_parent<<"\n";
  // std::cerr << "nodes="<<nodes<<"\n";
  nodes_of_cells.set_size(mesh.getNNodesPerCell(), cells_on_bdry.n_cols);
  bmi_cells_to_parent.set_size(2,cells_on_bdry.n_cols);
  for(int i = 0; i < cells_on_bdry.n_cols; i++)
  {
    int iK = cells_on_bdry(0,i);
    int iS = cells_on_bdry(1,i);
    int iil = cells_on_bdry(2,i);
    for(int ii=0;ii<n_nodes_per_side_parent;ii++)
    {
      int iN = nodes_of_sides_parents(ii,iS);
      nodes_of_cells(ii, i) = nodes_global_to_local[iN];
    }
    bmi_cells_to_parent(0,i) = iK;
    bmi_cells_to_parent(1,i) = iil;
  }
  // std::cerr << "nodes_of_cells="<<nodes_of_cells<<"\n";
  // std::cerr << "color_to_bdry2="<<color_to_bdry2<<"\n";

  alat::IntMap sizesOfColor;
  for(alat::Map<int, alat::armaimat>::const_iterator p=color_to_bdry2.begin();p!=color_to_bdry2.end();p++)
  {
    for(int i=0;i<p->second.n_cols;i++)
    {
      bool present = true;
      for(int ii=0;ii<p->second.n_rows;ii++)
      {
        if(nodes_global_to_local.find(p->second(ii,i))==nodes_global_to_local.end())
        {
          present=false;
          break;
        }
      }
      if(present)
      {
        if(sizesOfColor.hasKey(p->first)) sizesOfColor[p->first]++;
        else sizesOfColor[p->first]=1;
      }
    }
  }
  // std::cerr << "sizesOfColor="<<sizesOfColor<<"\n";

  alat::Map<int, alat::armaimat> color_to_bdry_new;
  for(alat::Map<int, alat::armaimat>::const_iterator p=color_to_bdry2.begin();p!=color_to_bdry2.end();p++)
  {
    if(sizesOfColor.find(p->first)==sizesOfColor.end())
    {
      continue;
    }
    color_to_bdry_new[p->first].set_size(p->second.n_rows,sizesOfColor[p->first]);
    // std::cerr << "n_nodes_per_cell " << n_nodes_per_cell << "\n";
    // std::cerr << "p->second.n_rows " << p->second.n_rows << "\n";
    // assert(n_nodes_per_cell==p->second.n_rows);
    int count=0;
    for(int i=0;i<p->second.n_cols;i++)
    {
      bool present = true;
      for(int ii=0;ii<p->second.n_rows;ii++)
      {
        if(nodes_global_to_local.find(p->second(ii,i))==nodes_global_to_local.end())
        {
          present=false;
          break;
        }
      }
      if(present)
      {
        for(int ii=0;ii<p->second.n_rows;ii++)
        {
          color_to_bdry_new[p->first](ii,count) = nodes_global_to_local[p->second(ii,i)];
        }
        count++;
      }
    }
    // std::cerr << "?? what ?? " << p->first << " -> " << p->second << "\n";
  }
  // std::cerr << "color_to_bdry_new="<<color_to_bdry_new<<"\n";
  if(mesh.getDimension())
  {
    _callSidesConstructor(&mesh, color_to_bdry_new);
    // std::cerr << "BoundaryInformation = " << _boundaryinformationmap << "\n";
  }
  // else
  // {
  //   // assert(n_dim==0);
  // }
}

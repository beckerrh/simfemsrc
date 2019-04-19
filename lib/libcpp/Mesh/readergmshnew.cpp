#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Mesh/readergmshnew.hpp"
#include  <fstream>

using namespace mesh;

/*---------------------------------------------------------*/
ReaderGmshNew::~ReaderGmshNew() {}
ReaderGmshNew::ReaderGmshNew() : alat::InterfaceBase()
{
  _size_of_elem.set_size(16);
  _size_of_elem[1] = 2;
  _size_of_elem[2] = 3;
  _size_of_elem[3] = 4;
  _size_of_elem[4] = 4;
  _size_of_elem[5] = 8;
  _size_of_elem[6] = 6;
  _size_of_elem[7] = 5;
  _size_of_elem[8] = 3;
  _size_of_elem[9] = 6;
  _size_of_elem[10] = 9;
  _size_of_elem[11] = 10;
  _size_of_elem[12] = 27;
  _size_of_elem[13] = 18;
  _size_of_elem[14] = 15;
  _size_of_elem[15] = 1;
}

std::string ReaderGmshNew::getClassName() const
{
  return "ReaderGmshNew";
}

/*---------------------------------------------------------*/
void ReaderGmshNew::read(std::string filename)
{
  // std::cerr << "MeshBase::readGmsh() " <<filename << "\n";
  std::ifstream file(filename.c_str() );
  if(not file.is_open() )
  {
    _error_string("read", "cannot open file", filename);
  }

  std::string str;
  int format = 0, size = 0;
  double version;
  _binary = false;

  std::getline(file, str);
  if(not file)
  {
    _error_string("read", "cannot read file", filename);
  }
  if(str.find("$MeshFormat") != static_cast<std::string::size_type>( 0 ) )
  {
    _error_string("read", "cannot read ", str);
  }
  file >> version >> format >> size;
  std::getline(file, str);
  std::cerr << " version " << version<< " format " << format<< " size " << size << "\n";
  std::getline(file, str);
  std::cerr << "str " << str << "\n";
  if(version != 4.1)
  {
    _error_string("read", "unknown version", version);
  }
  if(format == 1)
  {
    _binary = true;
  }
  else if(format)
  {
    _error_string("read", "unknown format", format);
  }
  if(_binary)
  {
    // int one;
    // file.read( reinterpret_cast<char*>(&one), sizeof( int ) );
    // std::cerr << "one="<<one << "\n";
    // assert(one==1);
    std::getline(file, str);
  }
  while(true)
  {
    std::getline(file, str);
    std::cerr << "reading " << str << "\n";
    if(file.eof() )
    {
      break;
    }
    if(str.find("$Nodes") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_nodes_binaray(file);
      }
      else
      {
        _read_nodes_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else if(str.find("$Elements") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_elements_binaray(file);
      }
      else
      {
        _read_elements_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else if(str.find("$PhysicalNames") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_physicalnames_binaray(file);
      }
      else
      {
        _read_physicalnames_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else if(str.find("$Entities") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_entities_binaray(file);
      }
      else
      {
        _read_entities_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else if(str.find("$PartitionedEntities") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_partitionedentities_binaray(file);
      }
      else
      {
        _read_partitionedentities_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else if(str.find("$Periodic") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_periodic_binaray(file);
      }
      else
      {
        _read_periodic_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else if(str.find("$GhostElements") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_ghostelements_binaray(file);
      }
      else
      {
        _read_ghostelements_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else if(str.find("$NodeData") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_nodedata_binaray(file);
      }
      else
      {
        _read_nodedata_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else if(str.find("$ElementData") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_elementdata_binaray(file);
      }
      else
      {
        _read_elementdata_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else if(str.find("$ElementNodeData") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_elementnodedata_binaray(file);
      }
      else
      {
        _read_elementnodedata_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else if(str.find("$InterpolationScheme") == static_cast<std::string::size_type>( 0 ) )
    {
      if(_binary)
      {
        _read_interpolationscheme_binaray(file);
      }
      else
      {
        _read_interpolationscheme_ascii(file);
      }
      std::getline(file, str);
      std::getline(file, str);
    }
    else
    {
      _error_string("read", "unknown field", str);
    }
    continue;
  }
  file.close();
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_entities_ascii(std::ifstream& file)
{
  int numPoints, numCurves, numSurfaces, numVolumes;
  int nphystags;

  file >> numPoints >> numCurves >> numSurfaces >> numVolumes;
  _pointtag.tag.set_size(numPoints);
  _pointtag.coord.set_size(3, numPoints);
  _pointtag.physicaltag.set_size(numPoints);

  _curvetag.tag.set_size(numCurves);
  _curvetag.coordmin.set_size(3, numCurves);
  _curvetag.coordmax.set_size(3, numCurves);
  _curvetag.physicaltag.set_size(numCurves);
  _curvetag.physicaltaglower.set_size(numCurves);

  _surfacetag.tag.set_size(numSurfaces);
  _surfacetag.coordmin.set_size(3, numSurfaces);
  _surfacetag.coordmax.set_size(3, numSurfaces);
  _surfacetag.physicaltag.set_size(numSurfaces);
  _surfacetag.physicaltaglower.set_size(numSurfaces);

  _volumetag.tag.set_size(numVolumes);
  _volumetag.coordmin.set_size(3, numVolumes);
  _volumetag.coordmax.set_size(3, numVolumes);
  _volumetag.physicaltag.set_size(numVolumes);
  _volumetag.physicaltaglower.set_size(numVolumes);

  for(int i = 0; i < numPoints; i++)
  {
    file >> _pointtag.tag[i] >> _pointtag.coord(0, i) >> _pointtag.coord(1, i) >> _pointtag.coord(2, i) >> nphystags;
    _pointtag.physicaltag[i].set_size(nphystags);
    for(int ii = 0; ii < nphystags; ii++)
    {
      file >> _pointtag.physicaltag[i][ii];
    }
  }

  for(int i = 0; i < numCurves; i++)
  {
    file >> _curvetag.tag[i] >> _curvetag.coordmin(0, i) >> _curvetag.coordmin(1, i) >> _curvetag.coordmin(2, i) >> _curvetag.coordmax(0, i) >> _curvetag.coordmax(1, i) >> _curvetag.coordmax(2, i) >> nphystags;
    _curvetag.physicaltag[i].set_size(nphystags);
    for(int ii = 0; ii < nphystags; ii++)
    {
      file >> _curvetag.physicaltag[i][ii];
    }
    file >> nphystags;
    _curvetag.physicaltaglower[i].set_size(nphystags);
    for(int ii = 0; ii < nphystags; ii++)
    {
      file >> _curvetag.physicaltaglower[i][ii];
    }
  }

  for(int i = 0; i < numSurfaces; i++)
  {
    file >> _surfacetag.tag[i] >> _surfacetag.coordmin(0, i) >> _surfacetag.coordmin(1, i) >> _surfacetag.coordmin(2, i) >> _surfacetag.coordmax(0, i) >> _surfacetag.coordmax(1, i) >> _surfacetag.coordmax(2, i) >> nphystags;
    _surfacetag.physicaltag[i].set_size(nphystags);
    for(int ii = 0; ii < nphystags; ii++)
    {
      file >> _surfacetag.physicaltag[i][ii];
    }
    file >> nphystags;
    _surfacetag.physicaltaglower[i].set_size(nphystags);
    for(int ii = 0; ii < nphystags; ii++)
    {
      file >> _surfacetag.physicaltaglower[i][ii];
    }
  }

  for(int i = 0; i < numVolumes; i++)
  {
    file >> _volumetag.tag[i] >> _volumetag.coordmin(0, i) >> _volumetag.coordmin(1, i) >> _volumetag.coordmin(2, i) >> _volumetag.coordmax(0, i) >> _volumetag.coordmax(1, i) >> _volumetag.coordmax(2, i) >> nphystags;
    _volumetag.physicaltag[i].set_size(nphystags);
    for(int ii = 0; ii < nphystags; ii++)
    {
      file >> _volumetag.physicaltag[i][ii];
    }
    file >> nphystags;
    _volumetag.physicaltaglower[i].set_size(nphystags);
    for(int ii = 0; ii < nphystags; ii++)
    {
      file >> _volumetag.physicaltaglower[i][ii];
    }
  }
  // std::cerr << "_pointtag" << _pointtag.tag << "\n" << _pointtag.coord << "\n" << _pointtag.physicaltag << "\n";
// std::cerr << "_curvetag" << _curvetag.tag << "\n" << _curvetag.coordmin << "\n" << _curvetag.coordmax << "\n" << _curvetag.physicaltag << "\n" << _curvetag.physicaltaglower << "\n";
  // std::cerr << "_surfacetag" << _surfacetag.tag << "\n" << _surfacetag.coordmin << "\n" << _surfacetag.coordmax << "\n" << _surfacetag.physicaltag << "\n" << _surfacetag.physicaltaglower << "\n";
  // std::cerr << "_volumetag" << _volumetag.tag << "\n" << _volumetag.coordmin << "\n" << _volumetag.coordmax << "\n" << _volumetag.physicaltag << "\n" << _volumetag.physicaltaglower << "\n";
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_nodes_binaray(std::ifstream& file)
{
  _notWritten("_read_nodes_binaray");
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_nodes_ascii(std::ifstream& file)
{
  int numEntityBlocks, numNodes, minNodeTag, maxNodeTag;
  file >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;
  std::cerr << "numEntityBlocks " << numEntityBlocks << " numNodes " << numNodes << " minNodeTag " << minNodeTag << " maxNodeTag " << maxNodeTag << "\n";
  _nodes.set_size(3, numNodes);
  int entityDim, entityTag, parametric, numNodesBlock;
  int count = 0;
  for(int i = 0; i < numEntityBlocks; i++)
  {
    file >> entityDim >> entityTag >> parametric >> numNodesBlock;
    assert(parametric == 0);
// std::cerr << " entityDim " << entityDim << " ";
// std::cerr << " entityTag " << entityTag << " ";
// std::cerr << " parametric " << parametric << " ";
// std::cerr << " numNodesBlock " << numNodesBlock << " ";
    alat::armaivec nodeTag(numNodesBlock);
    for(int ii = 0; ii < numNodesBlock; ii++)
    {
      file >> nodeTag[ii];
    }
    // std::cerr << " nodeTag ";
    for(int ii = 0; ii < numNodesBlock; ii++)
    {
      // std::cerr << nodeTag[ii] << "\n";
      _nodeid2id[nodeTag[ii]] = count;
      file >> _nodes(0, count) >> _nodes(1, count) >> _nodes(2, count);
      count++;
    }
    // std::cerr << "\n";
  }
  // std::cerr << " _nodes " << _nodes.t() << " ";
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_elements_binaray(std::ifstream& file)
{
  _notWritten("_read_nodes_binaray");
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_elements_ascii(std::ifstream& file)
{
  alat::armaivec& col = _elems.col();
  alat::armaivec& rowstart = _elems.rowstart();

  int numEntityBlocks, numElements, minElementTag, maxElementTag;
  file >> numEntityBlocks >> numElements >> minElementTag >> maxElementTag;
  std::cerr << "numEntityBlocks numElements " << numEntityBlocks << " " << numElements << " "<< minElementTag << " " << maxElementTag<<"\n";
  int entityDim, entityTag, elementType, numElementsBlock;
  rowstart.set_size(numElements+1);
  rowstart[0] = 0;
  int count = 0;
  std::streampos elemstart = file.tellg();
  for(int i = 0; i < numEntityBlocks; i++)
  {
    file >> entityDim >> entityTag >> elementType >> numElementsBlock;
    std::cerr << "entityDim entityTag elementType numElementsBlock " << entityDim << " " << entityTag << " "<< elementType << " " << numElementsBlock<<"\n";
    int nnodes = _size_of_elem[elementType];
    alat::armaivec data(nnodes+1);
    for(int ii = 0; ii < numElementsBlock; ii++)
    {
      for(int iii = 0; iii < nnodes+1; iii++)
      {
        file >> data[iii];
      }
      int id = data[0];
      _cellid2id[id] = count;
      rowstart[count+1] = rowstart[count] + nnodes +2;
      count++;
    }
  }
  int colsize = rowstart[numElements];
  col.set_size(colsize);
  file.seekg(elemstart);
  colsize = 0;
  for(int i = 0; i < numEntityBlocks; i++)
  {
    file >> entityDim >> entityTag >> elementType >> numElementsBlock;
    int nnodes = _size_of_elem[elementType];
    alat::armaivec data(nnodes+1);
    for(int ii = 0; ii < numElementsBlock; ii++)
    {
      col[colsize] = elementType;
      for(int iii = 0; iii < nnodes+1; iii++)
      {
        file >> col[colsize+iii+1];
      }
      colsize += nnodes+2;
    }
  }
  // std::cerr << "_elems = " << _elems << "\n";
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_physicalnames_binaray(std::ifstream& file)
{
  _notWritten("_read_physicalnames_binaray");
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_physicalnames_ascii(std::ifstream& file)
{
  _notWritten("_read_physicalnames_ascii");
}

void ReaderGmshNew::_read_entities_binaray(std::ifstream& file)
{
  _notWritten("_read_entities_binaray");
}

void ReaderGmshNew::_read_partitionedentities_binaray(std::ifstream& file)
{
  _notWritten("_read_partitionedentities_binaray");
}

void ReaderGmshNew::_read_partitionedentities_ascii(std::ifstream& file)
{
  _notWritten("_read_partitionedentities_ascii");
}

void ReaderGmshNew::_read_periodic_binaray(std::ifstream& file)
{
  _notWritten("_read_periodic_binaray");
}

void ReaderGmshNew::_read_periodic_ascii(std::ifstream& file)
{
  _notWritten("_read_periodic_ascii");
}

void ReaderGmshNew::_read_ghostelements_binaray(std::ifstream& file)
{
  _notWritten("_read_ghostelements_binaray");
}

void ReaderGmshNew::_read_ghostelements_ascii(std::ifstream& file)
{
  _notWritten("_read_ghostelements_ascii");
}

void ReaderGmshNew::_read_nodedata_binaray(std::ifstream& file)
{
  _notWritten("_read_nodedata_binaray");
}

void ReaderGmshNew::_read_nodedata_ascii(std::ifstream& file)
{
  _notWritten("_read_nodedata_ascii");
}

void ReaderGmshNew::_read_elementdata_binaray(std::ifstream& file)
{
  _notWritten("_read_elementdata_binaray");
}

void ReaderGmshNew::_read_elementdata_ascii(std::ifstream& file)
{
  _notWritten("_read_elementdata_ascii");
}

void ReaderGmshNew::_read_elementnodedata_binaray(std::ifstream& file)
{
  _notWritten("_read_elementnodedata_binaray");
}

void ReaderGmshNew::_read_elementnodedata_ascii(std::ifstream& file)
{
  _notWritten("_read_elementnodedata_ascii");
}

void ReaderGmshNew::_read_interpolationscheme_binaray(std::ifstream& file)
{
  _notWritten("_read_interpolationscheme_binaray");
}

void ReaderGmshNew::_read_interpolationscheme_ascii(std::ifstream& file)
{
  _notWritten("_read_interpolationscheme_ascii");
}

/*---------------------------------------------------------*/
void ReaderGmshNew::setMesh(MeshUnitInterface* mesh)
{
  // std::cerr << "ReaderGmshNew::setMesh() _npartsmax="<<_npartsmax<<"\n";
// if(_npartsmax)
// {
//   _setMeshAndInterfaces(mesh);
//   return;
// }
  arma::mat& nodes = mesh->getNodesAndNodesOfCells()._nodes;
  alat::armaimat& nodes_of_cells = mesh->getNodesAndNodesOfCells()._nodes_of_cells;
  int gmshcell = mesh->getGmshCellType();
  int gmshside = mesh->getGmshSideType();
  int gmshsideside = mesh->getGmshSideSideType();
  int gmshsidesideside = mesh->getGmshSideSideSideType();
  int n_nodes_per_cell = mesh->getNNodesPerCell();
  int n_nodes_per_side = mesh->getNNodesPerSide();
  int n_dim = mesh->getDimension();
  int n_nodes_per_side_side = 1;
  if(n_dim == 3)
  {
    n_nodes_per_side_side = 2;
  }

  // std::cerr << "_elems " << _elems << "\n";
  // std::cerr << "_nodeid2id " << _nodeid2id << "\n";

  nodes = _nodes;

  // alat::Map<int, alat::armaimat> color_to_bdry1, color_to_bdry2, color_to_bdry3;

  if(n_dim == 3)
  {
    _volumetag.tag;
  }
  else
  {
    std::cerr << "_surfacetag.tag" << _surfacetag.tag << "\n";
    std::cerr << "_surfacetag.physicaltag" << _surfacetag.physicaltag << "\n";
    std::cerr << "_surfacetag.physicaltaglower" << _surfacetag.physicaltaglower << "\n";

    std::cerr << "_curvetag.tag" << _curvetag.tag << "\n";
    std::cerr << "_curvetag.physicaltag" << _curvetag.physicaltag << "\n";
    std::cerr << "_curvetag.physicaltaglower" << _curvetag.physicaltaglower << "\n";
  }
  assert(0);


  // First round : find sizes
  int countcells = 0;
  for(int i = 0; i < _elems.n(); i++)
  {
    int pos = _elems.rowstart(i);
    int type = _elems.col(pos);
    int id = _elems.col(pos+1);
    if(type == gmshcell)
    {
      countcells++;
    }
  }
  nodes_of_cells.set_size(n_nodes_per_cell, countcells);
  // Second round : save elements
  countcells = 0;
  for(int i = 0; i < _elems.n(); i++)
  {
    int pos = _elems.rowstart(i);
    int type = _elems.col(pos);
    int id = _elems.col(pos+1);
    if(type == gmshcell)
    {
      assert( _size_of_elem[type] == n_nodes_per_cell);
      for(int ii = 0; ii < n_nodes_per_cell; ii++)
      {
        nodes_of_cells(ii, countcells) = _nodeid2id[_elems.col(pos+2+ii)];
      }
      countcells++;
    }
  }
  // std::cerr << "nodes_of_cells " << nodes_of_cells.t()<<"\n";
  // std::cerr << "color_to_bdry1 " << color_to_bdry1<<"\n";
  // std::cerr << "color_to_bdry2 " << color_to_bdry2<<"\n";
  // std::cerr << "color_to_bdry3 " << color_to_bdry3<<"\n";
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_setMeshAndInterfaces(MeshUnitInterface* mesh)
{
  assert(0);
}

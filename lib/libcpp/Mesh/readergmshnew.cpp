#include  "Mesh/readergmshnew.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
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
std::string ReaderGmshNew::getClassName() const {return "ReaderGmshNew";}

/*---------------------------------------------------------*/
void ReaderGmshNew::read(std::string filename)
{
  // std::cerr << "MeshBase::readGmsh() " <<filename << "\n";
  std::ifstream file(filename.c_str());
  if(not file.is_open()) {_error_string("read","cannot open file", filename);}

  std::string str;
  int format=0, size=0;
  double version;
  _binary=false;

  std::getline(file, str);
  if (not file)
  {
    _error_string("read", "cannot read file",filename);
  }
  if (str.find("$MeshFormat") != static_cast<std::string::size_type>(0))
  {
    _error_string("read", "cannot read ",str);
  }
  file >> version >> format >> size;
  std::getline(file, str);
  // std::cerr << " version " << version<< " format " << format<< " size " << size << "\n";
  std::getline(file, str);
  // std::cerr << "str " << str << "\n";
  if(version != 2.2)
  {
    _error_string("read", "unknown version",version);
  }
  if (format==1)
  {
    _binary = true;
  }
  else if (format)
  {
    _error_string("read", "unknown format",format);
  }
  if(_binary)
  {
    // int one;
    // file.read( reinterpret_cast<char*>(&one), sizeof( int ) );
    // std::cerr << "one="<<one << "\n";
    // assert(one==1);
    std::getline(file, str);
  }
  while (true)
  {
    std::getline(file, str);
    // std::cerr << "reading " << str << "\n";
    if (file.eof())
    {
      break;
    }
    if (str.find("$Nodes") == static_cast<std::string::size_type>(0))
    {
      if(_binary) {_read_nodes_binaray(file);}
      else {_read_nodes_ascii(file);}
      std::getline(file, str);
      std::getline(file, str);
    }
    else if (str.find("$Elements") == static_cast<std::string::size_type>(0))
    {
      if(_binary) {_read_elements_binaray(file);}
      else {_read_elements_ascii(file);}
      std::getline(file, str);
      std::getline(file, str);
    }
    else
    {
      _error_string("read", "unknown field",str);
    }
    continue;
  }
  file.close();
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_nodes_binaray(std::ifstream& file)
{
  std::string str;
  int num_nodes, id;
  file >> num_nodes;
  std::getline(file, str);
  _nodes.set_size(3,num_nodes);
  // char buffer[100];
  for (int i=0; i<num_nodes; ++i)
  {
    file.read( reinterpret_cast<char*>(&id), sizeof( int ) );
    file.read( reinterpret_cast<char*>(&_nodes(0,i)), sizeof( double ) );
    file.read( reinterpret_cast<char*>(&_nodes(1,i)), sizeof( double ) );
    file.read( reinterpret_cast<char*>(&_nodes(2,i)), sizeof( double ) );
    _nodeid2id[id] = i;
  }
  // std::cerr << "_nodes="<<_nodes.t()<<"\n";
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_nodes_ascii(std::ifstream& file)
{
  std::string str;
  int num_nodes, id;
  file >> num_nodes;
  std::getline(file, str);
  _nodes.set_size(3,num_nodes);
  double x, y, z;
  for (int i=0; i<num_nodes; ++i)
  {
    file >> id >> _nodes(0,i) >> _nodes(1,i) >> _nodes(2,i);
    _nodeid2id[id] = i;
  }
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_elements_binaray(std::ifstream& file)
{
  alat::armaivec& col = _elems.col();
  alat::armaivec& rowstart = _elems.rowstart();

  std::string str;
  int nelems;
  file >> nelems;
  std::getline(file, str);

  std::streampos elemstart=file.tellg();

  int id, type, ntags, physical, elementary, nparts, num_elm_follow;
  int buffer[10];
  _npartsmax=0;
  alat::armaivec rowsizes(nelems);
  int header[3];
  for(int i=0; i<nelems; i++)
  {
    // file >> id >> type >> ntags >> physical >> elementary;
    file.read( reinterpret_cast<char*>(&header), 3*sizeof( int ) );
    type = header[0];
    num_elm_follow = header[1];
    ntags = header[2];
    int ndata = 1+ntags+_size_of_elem[type];
    int data[ndata];
    file.read( reinterpret_cast<char*>(&data), ndata*sizeof( int ) );
    id = data[0];
    physical = data[1];
    elementary = data[2];
    // rowsizes[i] = ndata + 3;
    rowsizes[i] = ndata + 2;
    // std::cerr << "header = ";for(int ii=0;ii<3;ii++) std::cerr<< header[ii]<<" "; std::cerr << "\t";
    // std::cerr << "data = ";for(int ii=0;ii<ndata;ii++) std::cerr << data[ii]<< " "; std::cerr << "\n";
    // std::cerr << "id >> type >> ntags >> physical >> elementary " << id << " " << type << " "<< ntags<< " " << physical<< " " << elementary<<"\n";
    // std::cerr << "num_elm_follow = " << num_elm_follow << "\n";
    if(ntags>2)
    {
      nparts = data[3];
    }
    else
    {
      nparts=0;
    }
    // std::cerr << "nparts="<<nparts<< " ntags="<<ntags<<"\n";
    _npartsmax = std::max(_npartsmax,nparts);

    _cellid2id[id]=i;
  }
  _elems.set_size(nelems, arma::sum(rowsizes));
  file.seekg(elemstart);
  rowstart[0] = 0;
  int count=0;
  for(int i=0; i<nelems; ++i)
  {
    rowstart[i+1] = rowstart[i] + rowsizes[i];
    file.read( reinterpret_cast<char*>(&header), 3*sizeof( int ) );
    type = header[0];
    ntags = header[2];
    int ndata = 1+ntags+_size_of_elem[type];
    int data[ndata];
    file.read( reinterpret_cast<char*>(&data), ndata*sizeof( int ) );
    col[count++] = data[0];
    col[count++] = header[0];
    col[count++] = header[2];
    for(int ii=1;ii<ndata;ii++)  {col[count++] = data[ii];}
  }
}

/*---------------------------------------------------------*/
void ReaderGmshNew::_read_elements_ascii(std::ifstream& file)
{
  alat::armaivec& col = _elems.col();
  alat::armaivec& rowstart = _elems.rowstart();

  std::string str;
  int nelems;
  file >> nelems;
  std::getline(file, str);

  std::streampos elemstart=file.tellg();

  // std::cerr << "nelems " << nelems <<"\n";

  int id, type, ntags, physical, elementary, nparts;
  _npartsmax=0;
  alat::armaivec rowsizes(nelems);
  alat::armaivec header(3);
  // int header[3];
  for(int i=0; i<nelems; i++)
  {
    // file >> id >> type >> ntags >> physical >> elementary;
    for(int ii=0;ii<3;ii++) file >> header[ii];
    // header.load(file);
    id = header[0];
    type = header[1];
    ntags = header[2];

    // std::cerr << "id >> type >> ntags " << id << " " << type << " " << ntags <<"\n";
    // assert(ntags==2);

    int ndata = ntags+_size_of_elem[type];
    // std::cerr << "ndata " << ndata <<"\n";
    // int data[ndata];
    alat::armaivec data(ndata);
    for(int ii=0;ii<ndata;ii++) file >> data[ii];
    // file >> data;
    // data.load(file);
    // std::cerr << "id >> type >> ntags >> data " << id << " " << type << " " << ntags << " " << data.t() <<"\n";
    rowsizes[i] = ndata + 3;

    if(ntags>2)
    {
      nparts = data[2];
    }
    else
    {
      nparts=0;
    }
    _npartsmax = std::max(_npartsmax,nparts);
    _cellid2id[id]=i;
  }
  _elems.set_size(nelems, arma::sum(rowsizes));
  file.seekg(elemstart);
  rowstart[0] = 0;
  int count=0;
  for(int i=0; i<nelems; ++i)
  {
    rowstart[i+1] = rowstart[i] + rowsizes[i];
    for(int ii=0;ii<rowsizes[i];ii++)
    {
      file >> col[count++];
    }
  }
}

/*---------------------------------------------------------*/
void ReaderGmshNew::setMesh(MeshUnitInterface* mesh)
{
  // std::cerr << "ReaderGmshNew::setMesh() _npartsmax="<<_npartsmax<<"\n";
  if(_npartsmax)
  {
    _setMeshAndInterfaces(mesh);
    return;
  }
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
  if(n_dim==3) {n_nodes_per_side_side = 2;}

  // std::cerr << "_elems " << _elems << "\n";
  // std::cerr << "_nodeid2id " << _nodeid2id << "\n";

  nodes = _nodes;

  // First round : find sizes
  int countcells=0;
  alat::IntMap color_to_sizebdry1, color_to_sizebdry2, color_to_sizebdry3;
  for(int i=0;i<_elems.n(); i++)
  {
    int pos = _elems.rowstart(i);
    int type = _elems.col(pos+1);
    int physical = _elems.col(pos+3);
    if(type==gmshcell) {countcells++;}
    else if(type==gmshside)
    {
      if(color_to_sizebdry1.hasKey(physical)) color_to_sizebdry1[physical]++;
      else color_to_sizebdry1[physical]=1;
    }
    else if(type==gmshsideside)
    {
      if(color_to_sizebdry2.hasKey(physical)) color_to_sizebdry2[physical]++;
      else color_to_sizebdry2[physical]=1;
    }
    else if(type==gmshsidesideside)
    {
      if(color_to_sizebdry3.hasKey(physical)) color_to_sizebdry3[physical]++;
      else color_to_sizebdry3[physical]=1;
    }
  }
  nodes_of_cells.set_size(n_nodes_per_cell, countcells);
  alat::IntMap bdry1count, bdry2count, bdry3count;
  alat::IntMap::const_iterator p;
  for(p=color_to_sizebdry1.begin(); p!=color_to_sizebdry1.end();p++)
  {
    color_to_bdry1[p->first].set_size(n_nodes_per_side, p->second);
    bdry1count[p->first]=0;
  }
  for(p=color_to_sizebdry2.begin(); p!=color_to_sizebdry2.end();p++)
  {
    color_to_bdry2[p->first].set_size(n_nodes_per_side_side, p->second);
    bdry2count[p->first]=0;
  }
  for(p=color_to_sizebdry3.begin(); p!=color_to_sizebdry3.end();p++)
  {
    color_to_bdry3[p->first].set_size(1, p->second);
    bdry3count[p->first]=0;
  }

  // Second round : save elements
  countcells=0;
  for(int i=0;i<_elems.n(); i++)
  {
    int pos = _elems.rowstart(i);
    int id = _elems.col(pos);
    int type = _elems.col(pos+1);
    int ntags = _elems.col(pos+2);
    assert(ntags==2);
    int physical = _elems.col(pos+3);
    int elementary = _elems.col(pos+4);
    if(type==gmshcell)
    {
      assert( _size_of_elem[type] == n_nodes_per_cell);
      for(int ii=0;ii<n_nodes_per_cell;ii++)
      {
        nodes_of_cells(ii,countcells) = _nodeid2id[_elems.col(pos+5+ii)];
      }
      countcells++;
    }
    else if(type==gmshside)
    {
      assert( _size_of_elem[type] == n_nodes_per_side);
      for(int ii=0;ii<n_nodes_per_side;ii++)
      {
        color_to_bdry1[physical](ii,bdry1count[physical])=_nodeid2id[_elems.col(pos+5+ii)];
      }
      bdry1count[physical]++;
    }
    else if(type==gmshsideside)
    {
      assert( _size_of_elem[type] == n_nodes_per_side_side);
      for(int ii=0;ii<n_nodes_per_side_side;ii++)
      {
        color_to_bdry2[physical](ii,bdry2count[physical])=_nodeid2id[_elems.col(pos+5+ii)];
      }
      bdry2count[physical]++;
    }
    else if(type==gmshsidesideside)
    {
      assert( _size_of_elem[type] == 1);
      for(int ii=0;ii<1;ii++)
      {
        color_to_bdry3[physical](ii,bdry3count[physical])=_nodeid2id[_elems.col(pos+5+ii)];
      }
      bdry3count[physical]++;
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
  if(n_dim==3) {n_nodes_per_side_side = 2;}

  nodes = _nodes;

  // First round : find sizes
  int countcells=0;
  alat::IntMap color_to_sizebdry1, color_to_sizebdry2, color_to_sizebdry3;
  alat::IntMap partion_to_sizecells;
  for(int i=0;i<_elems.n(); i++)
  {
    // std::cerr << "_elem(i) = ";
    // for(int pos2=_elems.rowstart(i); pos2 < _elems.rowstop(i); pos2++)
    // {
    //   std::cerr << " " << _elems.col(pos2);
    // }
    // std::cerr << "\n";
    int pos = _elems.rowstart(i);
    int type = _elems.col(pos+1);
    int ntags = _elems.col(pos+2);
    int physical = _elems.col(pos+3);
    if(ntags<=2)
    {
      _error_string("_setMeshAndInterfaces", "unexpected ntags=2");
    }
    int nparts = _elems.col(pos+5);
    assert(ntags=nparts+3);
    alat::armaivec parts(nparts);
    for(int ii=0;ii<nparts;ii++) {parts[ii] = _elems.col(pos+6+ii);}
    // std::cerr << "i nparts ntags physical parts " << i << " " << nparts<< " " << ntags << " " << physical << " " << parts.t();
    if(type==gmshcell)
    {
      countcells++;
      for(int ii=0;ii<nparts;ii++)
      {
        if(parts[ii]<0)
        {
          int partition = parts[ii];
          if(partion_to_sizecells.hasKey(partition)) partion_to_sizecells[partition]++;
          else partion_to_sizecells[partition]=1;
        }
      }
    }
    else if(type==gmshside)
    {
      if(color_to_sizebdry1.hasKey(physical)) color_to_sizebdry1[physical]++;
      else color_to_sizebdry1[physical]=1;
    }
    else if(type==gmshsideside)
    {
      if(color_to_sizebdry2.hasKey(physical)) color_to_sizebdry2[physical]++;
      else color_to_sizebdry2[physical]=1;
    }
    else if(type==gmshsidesideside)
    {
      if(color_to_sizebdry3.hasKey(physical)) color_to_sizebdry3[physical]++;
      else color_to_sizebdry3[physical]=1;
    }
  }
  // std::cerr << "countcells=" << countcells << "\n";
  // std::cerr << "color_to_sizebdry1=" << color_to_sizebdry1 << "\n";
  // std::cerr << "color_to_sizebdry2=" << color_to_sizebdry2 << "\n";
  // std::cerr << "color_to_sizebdry3=" << color_to_sizebdry3 << "\n";
  // std::cerr << "partion_to_sizecells=" << partion_to_sizecells << "\n";

  nodes_of_cells.set_size(n_nodes_per_cell, countcells);
  alat::IntMap bdry1count, bdry2count, bdry3count, partscount;
  alat::IntMap::const_iterator p;
  for(p=color_to_sizebdry1.begin(); p!=color_to_sizebdry1.end();p++)
  {
    color_to_bdry1[p->first].set_size(n_nodes_per_side, p->second);
    bdry1count[p->first]=0;
  }
  for(p=color_to_sizebdry2.begin(); p!=color_to_sizebdry2.end();p++)
  {
    color_to_bdry2[p->first].set_size(n_nodes_per_side_side, p->second);
    bdry2count[p->first]=0;
  }
  for(p=color_to_sizebdry3.begin(); p!=color_to_sizebdry3.end();p++)
  {
    color_to_bdry3[p->first].set_size(1, p->second);
    bdry3count[p->first]=0;
  }
  for(p=partion_to_sizecells.begin(); p!=partion_to_sizecells.end();p++)
  {
    partition_to_cells[p->first].set_size(p->second);
    partscount[p->first]=0;
  }

  // Second round : save elements
  countcells=0;
  for(int i=0;i<_elems.n(); i++)
  {
    int pos = _elems.rowstart(i);
    int id = _elems.col(pos);
    int type = _elems.col(pos+1);
    int ntags = _elems.col(pos+2);
    int physical = _elems.col(pos+3);
    int elementary = _elems.col(pos+4);
    int nparts = _elems.col(pos+5);
    assert(ntags=nparts+3);
    alat::armaivec parts(nparts);
    for(int ii=0;ii<nparts;ii++) {parts[ii] = _elems.col(pos+6+ii);}
    int poscell = pos + ntags + 3;

    if(type==gmshcell)
    {
      assert( _size_of_elem[type] == n_nodes_per_cell);
      for(int ii=0;ii<n_nodes_per_cell;ii++)
      {
        nodes_of_cells(ii,countcells) = _nodeid2id[_elems.col(poscell+ii)];
      }
      for(int ii=0;ii<nparts;ii++)
      {
        if(parts[ii]<0)
        {
          int partition = parts[ii];
          partition_to_cells[partition][partscount[partition]++] = countcells;
          if(partion_to_sizecells.hasKey(partition)) partion_to_sizecells[partition]++;
          else partion_to_sizecells[partition]=1;
        }
      }
      countcells++;
    }
    else if(type==gmshside)
    {
      assert( _size_of_elem[type] == n_nodes_per_side);
      for(int ii=0;ii<n_nodes_per_side;ii++)
      {
        color_to_bdry1[physical](ii,bdry1count[physical])=_nodeid2id[_elems.col(poscell+ii)];
      }
      bdry1count[physical]++;
    }
    else if(type==gmshsideside)
    {
      assert( _size_of_elem[type] == n_nodes_per_side_side);
      for(int ii=0;ii<n_nodes_per_side_side;ii++)
      {
        color_to_bdry2[physical](ii,bdry2count[physical])=_nodeid2id[_elems.col(poscell+ii)];
      }
      bdry2count[physical]++;
    }
    else if(type==gmshsidesideside)
    {
      assert( _size_of_elem[type] == 1);
      for(int ii=0;ii<1;ii++)
      {
        color_to_bdry3[physical](ii,bdry3count[physical])=_nodeid2id[_elems.col(poscell+ii)];
      }
      bdry3count[physical]++;
    }
  }
  // std::cerr << "color_to_bdry1=" << color_to_bdry1 << "\n";
  // std::cerr << "color_to_bdry2=" << color_to_bdry2 << "\n";
  // std::cerr << "color_to_bdry3=" << color_to_bdry3 << "\n";
  // std::cerr << "partition_to_cells=" << partition_to_cells << "\n";
  // std::cerr << "nodes_of_cells=" << nodes_of_cells << "\n";
}

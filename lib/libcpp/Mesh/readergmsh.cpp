#include  "Mesh/readergmsh.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"

using namespace mesh;

/*---------------------------------------------------------*/
ReaderGmsh::~ReaderGmsh() {}
ReaderGmsh::ReaderGmsh(MeshUnitInterface* mesh) : alat::InterfaceBase(), _mesh(mesh) {}
std::string ReaderGmsh::getClassName() const {return "ReaderGmsh";}

/*---------------------------------------------------------*/
void ReaderGmsh::read(std::string filename)
{
  arma::mat& nodes = _mesh->getNodesAndNodesOfCells()._nodes;
  alat::armaimat& nodes_of_cells = _mesh->getNodesAndNodesOfCells()._nodes_of_cells;
  int gmshcell = _mesh->getGmshCellType();
  int gmshside = _mesh->getGmshSideType();
  int gmshsideside = _mesh->getGmshSideSideType();
  int gmshsidesideside = _mesh->getGmshSideSideSideType();
  int n_nodes_per_cell = _mesh->getNNodesPerCell();
  int n_nodes_per_side = _mesh->getNNodesPerSide();
  int n_dim = _mesh->getDimension();
  int n_nodes_per_side_side = 1;
  if(n_dim==3) {n_nodes_per_side_side = 2;}

  // std::cerr << "MeshBase::readGmsh() " <<filename << "\n";
  std::ifstream file(filename.c_str());
  assert(file.is_open());

  std::map<unsigned int, unsigned int> nodeid2id, cellid2id;
  std::string str;
  int format=0, size=0;
  double version;
  bool binary=false;

  std::getline(file, str);
  if (not file)
  {
    _error_string("readGmsh", "cannot read file",filename);
  }
  if (str.find("$MeshFormat") != static_cast<std::string::size_type>(0))
  {
    _error_string("readGmsh", "cannot read ",str);
  }
  file >> version >> format >> size;
  // std::cerr << " version " << version<< " format " << format<< " size " << size << "\n";
  if(version != 2.2)
  {
    _error_string("readGmsh", "unknown version",version);
  }
  if (format==1)
  {
    binary = true;
  }
  else if (format)
  {
    _error_string("readGmsh", "unknown format",format);
  }
  while (true)
  {
    std::getline(file, str);
    // std::cerr << "reading " << str << "\n";
    if (file.eof())
    {
      break;
    }
    if (str.find("$PhysicalNames") == static_cast<std::string::size_type>(0))
    {
      assert(not binary);
      unsigned int num_physical_groups;
      file >> num_physical_groups;
      std::getline(file, str);
      for (unsigned int i=0; i<num_physical_groups; ++i)
      {
        std::getline(file, str);
      }
    }

    else if (str.find("$Nodes") == static_cast<std::string::size_type>(0))
    {
      unsigned int num_nodes = 0;
      file >> num_nodes;
      std::getline(file, str);
      nodes.set_size(3,num_nodes);
      double x, y, z;
      unsigned int id;
      // std::cerr << "num_nodes = " << num_nodes << "\n";
      char buffer[100];
      for (unsigned int i=0; i<num_nodes; ++i)
      {
        if(binary)
        {
          file.read( reinterpret_cast<char*>(&id), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&nodes(0,i)), sizeof( double ) );
          file.read( reinterpret_cast<char*>(&nodes(1,i)), sizeof( double ) );
          file.read( reinterpret_cast<char*>(&nodes(2,i)), sizeof( double ) );
        }
        else
        {
          file >> id >> x >> y >> z;
          nodes(0,i) = x;
          nodes(1,i) = y;
          nodes(2,i) = z;
        }
        // std::cerr << "id = " << id << " nodes[i] = " << nodes[i] << "\n";
        nodeid2id[id] = i;
      }
      // std::cerr << "nodes = " << nodes << "\n";
    }
    else if (str.find("$Elements") == static_cast<std::string::size_type>(0))
    {
      unsigned int ntotalcells = 0, node_id = 0;
      file >> ntotalcells;
      std::getline(file, str);
      // std::cerr << "ntotalcells = " << ntotalcells << "\n";

      std::streampos elemstart=file.tellg();
      unsigned int id, type, physical, elementary, ntags;
      int tag;
      int countcells=0, countsides=0, countsidesides=0, countsidesidesides=0;
      // alat::armaivec buffer(n_nodes_per_cell);
      int buffer[n_nodes_per_cell];
      int toto;
      int nparts;
      int parts[3];
      alat::Map<int, int> color_to_sizebdry1, color_to_sizebdry2, color_to_sizebdry3;
      for (unsigned int iK=0; iK<ntotalcells; ++iK)
      {
        if(binary)
        {
          file.read( reinterpret_cast<char*>(&type), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&toto), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&ntags), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&id), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&physical), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&elementary), sizeof( int ) );
          if(type==gmshcell)
          {
            file.read( reinterpret_cast<char*>(buffer), n_nodes_per_cell*sizeof( int ) );
          }
          else if(type==gmshside)
          {
            file.read( reinterpret_cast<char*>(buffer), n_nodes_per_side*sizeof( int ) );
          }
          else if(type==gmshsideside)
          {
            file.read( reinterpret_cast<char*>(buffer), n_nodes_per_side_side*sizeof( int ) );
          }
          else if(type==gmshsidesideside)
          {
            file.read( reinterpret_cast<char*>(buffer), 1*sizeof( int ) );
          }
          else
          {
            _error_string("readGmsh", "unexpected type =",type);
          }
         // std::cerr << " id = " << id<< " type = " << type<< " ntags = " << ntags<< " toto = " << toto << " physical = " << physical << " elementary = " << elementary << " buffer = " << *buffer << "\n";
        }
        else
        {
          file >> id >> type >> ntags >> physical >> elementary;
          if(ntags>2)
          {
            file >> nparts;
            for(int ii=0;ii<nparts;ii++) file >> parts[ii];
          }
          else
          {
            nparts=-1;
          }
          std::getline(file, str);
        }
        if (ntags!= 2)
        {
          std::cerr << "ntags = " << ntags << " nparts " << nparts << "\n";
          std::cerr << "id >> type >> ntags >> physical >> elementary " << id << " " << type << " "<< ntags<< " " << physical<< " " << elementary<<"\n";
          _error_string("readGmsh", "expect ntags=2, but ntags=",ntags);
        }
        if(type==gmshcell) countcells++;
        else if(type==gmshside)
        {
          countsides++;
          if(color_to_sizebdry1.hasKey(physical)) color_to_sizebdry1[physical]++;
          else color_to_sizebdry1[physical]=1;
        }
        else if(type==gmshsideside)
        {
          countsidesides++;
          if(color_to_sizebdry2.hasKey(physical)) color_to_sizebdry2[physical]++;
          else color_to_sizebdry2[physical]=1;
        }
        else if(type==gmshsidesideside)
        {
          countsidesidesides++;
          if(color_to_sizebdry3.hasKey(physical)) color_to_sizebdry3[physical]++;
          else color_to_sizebdry3[physical]=1;
        }
      }
      // std::cerr << "gmshcell, gmshside, gmshsideside, gmshsidesideside="<<gmshcell<<" "<<gmshside<<" "<<gmshsideside<<" " << gmshsidesideside<<"\n";
      // std::cerr << "countcells, countsides="<<countcells<<" "<<countsides<<"\n";
      // _cells.set_size(countcells);

      nodes_of_cells.set_size(n_nodes_per_cell, countcells);
      // std::cerr << "\n" << "color_to_sizebdry1="<<color_to_sizebdry1;
      // std::cerr << "\n" << "color_to_sizebdry2="<<color_to_sizebdry2;
      // std::cerr << "\n" << "color_to_sizebdry3="<<color_to_sizebdry3;
      alat::Map<int,int> bdry1count, bdry2count, bdry3count;
      alat::Map<int,int>::const_iterator p;
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

      file.seekg(elemstart);
      countcells=0;
      int nodeid;
      for (unsigned int iK=0; iK<ntotalcells; ++iK)
      {
        if(binary)
        {
          file.read( reinterpret_cast<char*>(&type), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&toto), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&ntags), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&id), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&physical), sizeof( int ) );
          file.read( reinterpret_cast<char*>(&elementary), sizeof( int ) );
          // std::cerr << " id = " << id<< " type = " << type<< " ntags = " << ntags<< " toto = " << toto << " physical = " << physical << " elementary = " << elementary << " buffer = " << *buffer << "\n";
        }
        else
        {
          file >> id >> type >> ntags >> physical >> elementary;
        }
        if(type==gmshcell)
        {
          cellid2id[id] = countcells;
          for (unsigned int j = 0; j < n_nodes_per_cell; j++)
          {
            if(binary)
            {
              file.read( reinterpret_cast<char*>(&nodeid), sizeof( int ) );
            }
            else
            {
              file >> nodeid;
            }
            // std::cerr << "nodeid="<<nodeid<< " " << nodeid2id[nodeid] << "\n";
            // _cells[countcells][j] = nodeid2id[nodeid];
            nodes_of_cells(j, countcells) = nodeid2id[nodeid];
          }
          countcells++;
        }
        else if(type==gmshside)
        {
         // std::cerr << "id, type, ntags, physical, elementary = " << id<<" "<<type<<" "<<ntags<<" "<<physical<<" "<<elementary << " --" << n_nodes_per_side<<"\n";
          for (unsigned int j = 0; j < n_nodes_per_side; j++)
          {
            if(binary)
            {
              file.read( reinterpret_cast<char*>(&nodeid), sizeof( int ) );
            }
            else
            {
              file >> nodeid;
            }
            color_to_bdry1[physical](j,bdry1count[physical])=nodeid2id[nodeid];
          }
          bdry1count[physical]++;
        }
        // std::getline(file, str);
        // std::cerr << "???????????????? " << s << "\n";
        else if(type==gmshsideside)
        {
         // std::cerr << "id, type, ntags, physical, elementary = " << id<<" "<<type<<" "<<ntags<<" "<<physical<<" "<<elementary << " --" << n_nodes_per_side<<"\n";
          for (unsigned int j = 0; j < n_nodes_per_side_side; j++)
          {
            if(binary)
            {
              file.read( reinterpret_cast<char*>(&nodeid), sizeof( int ) );
            }
            else
            {
              file >> nodeid;
            }
            color_to_bdry2[physical](j,bdry2count[physical])=nodeid2id[nodeid];
          }
          bdry2count[physical]++;
        }
        // std::getline(file, str);
        // std::cerr << "???????????????? " << s << "\n";
        else if(type==gmshsidesideside)
        {
         // std::cerr << "id, type, ntags, physical, elementary = " << id<<" "<<type<<" "<<ntags<<" "<<physical<<" "<<elementary << " --" << n_nodes_per_side<<"\n";
          for (unsigned int j = 0; j < 1; j++)
          {
            if(binary)
            {
              file.read( reinterpret_cast<char*>(&nodeid), sizeof( int ) );
            }
            else
            {
              file >> nodeid;
            }
            color_to_bdry3[physical](j,bdry3count[physical])=nodeid2id[nodeid];
          }
          bdry3count[physical]++;
        }
        // std::getline(file, str);
        // std::cerr << "???????????????? " << s << "\n";
      }
    }
    continue;
  }
  file.close();

  // std::cerr << "@ color_to_bdry1 " << color_to_bdry1 << "\n";
  // std::cerr << "@ color_to_bdry2 " << color_to_bdry2 << "\n";
  // std::cerr << "@ color_to_bdry3 " << color_to_bdry3 << "\n";

}

#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Mesh/edgesconstructor.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Mesh/sidesandcells.hpp"
#include  <algorithm>
#include  <limits>

using namespace mesh;

/*---------------------------------------------------------*/
EdgesConstructor::~EdgesConstructor(){}
EdgesConstructor::EdgesConstructor(MeshUnitInterface* mesh) : _mesh(mesh) {}

/*---------------------------------------------------------*/
void EdgesConstructor::constructEdgesFromCells(bool debug)
{
  const alat::armaimat&  nodes_of_cells = _mesh->getNodesAndNodesOfCells().getNodesOfCells();
  const alat::armaimat&  nodes_of_sides = _mesh->getSidesAndCells().getNodesOfSides();
  const alat::armaimat&  sides_of_cells = _mesh->getSidesAndCells().getSidesOfCells();
  const alat::armaimat&  cells_of_sides = _mesh->getSidesAndCells().getCellsOfSides();

  alat::armaimat& nodes_of_edges = _mesh->getEdgesAndCells().getNodesOfEdges();
  alat::armaimat& edges_of_cells = _mesh->getEdgesAndCells().getEdgesOfCells();
  alat::armaicube& localnodes_of_edges_in_cells =_mesh->getEdgesAndCells().getLocalnodesOfEdgesOfCells();

  int nnodespercell = _mesh->getNNodesPerCell();
  int nedgespercell = _mesh->getNEdgesPerCell();
  int nsidespercell = _mesh->getNSidesPerCell();
  int ncells =  nodes_of_cells.n_cols;

  alat::Map<FixEdge, int> contcells;
  for(int iK=0; iK < ncells; iK++)
  {
    if(debug)
    {
      std::cerr << iK << " ==> " << alat::armaivec(nodes_of_cells.col(iK)).t();
    }
    for(int ii=0; ii < nedgespercell; ii++)
    {
      FixEdge s(_mesh->_getEdgeOfCell(iK, ii));
      if(debug)
      {
        std::cerr << "s=" << s << "\n";
      }
      std::sort( s.begin(), s.end() );
      alat::Map<FixEdge, int>::iterator p = contcells.find(s);
      if( p == contcells.end() )
      {
        int count=1;
        contcells.insert( std::make_pair(s, count) );
      }
      else
      {
        p->second++;
      }
    }
  }
  if(debug)
  {
    std::cerr << "nedges=" << contcells.size()<< "\n";
    std::cerr << "contcells=" << contcells<< "\n";
  }

  EdgeToInfo edgeinfo;
  for(alat::Map<FixEdge, int>::iterator p=contcells.begin();p!=contcells.end();p++)
  {
    edgeinfo.insert(std::make_pair(p->first, alat::armaimat(2,p->second)));
    p->second=0;
  }
  for(int iK=0; iK < ncells; iK++)
  {
    for(int ii=0; ii < nedgespercell; ii++)
    {
      FixEdge s(_mesh->_getEdgeOfCell(iK, ii));
      std::sort( s.begin(), s.end() );
      alat::Map<FixEdge, int>::iterator p = contcells.find(s);
      assert( p != contcells.end() );
      assert( edgeinfo.find(s) != edgeinfo.end() );
      edgeinfo[s](0, p->second) = iK;
      edgeinfo[s](1, p->second) = ii;
      p->second++;
    }
  }
  // if(debug)
  // {
  //   std::cerr << "nedges=" << edgeinfo.size()<< "\n";
  //   std::cerr << "edgeinfo=" << edgeinfo<< "\n";
  // }
  int nedges = contcells.size();
  nodes_of_edges.set_size(2, nedges);
  nodes_of_edges.fill(-1);
  edges_of_cells.set_size(nedgespercell, ncells);
  edges_of_cells.fill(-1);
  localnodes_of_edges_in_cells.set_size(2, nedgespercell, ncells);
  localnodes_of_edges_in_cells.fill(-1);

  int count=0;
  for(EdgeToInfo::const_iterator p=edgeinfo.begin();p!=edgeinfo.end();p++)
  {
    // std::cerr << "count=" << count << " p->first=" << p->first << "\n";
    nodes_of_edges(0, count) = p->first[0];
    nodes_of_edges(1, count) = p->first[1];
    const alat::armaimat& edgeinfo = p->second;
    for(int ii=0;ii<edgeinfo.n_cols;ii++)
    {
      int iK = edgeinfo(0,ii);
      // std::cerr << "iK=" << iK << " nodes_of_cells="<< alat::armaivec(nodes_of_cells.col(iK)).t();
      int iiE = edgeinfo(1,ii);
      edges_of_cells(iiE, iK) = count;
      for(int iii=0;iii<nnodespercell;iii++)
      {
        if(p->first[0] == nodes_of_cells(iii,iK))
        {
          localnodes_of_edges_in_cells(0, iiE, iK) = iii;
        }
        else if(p->first[1] == nodes_of_cells(iii,iK))
        {
          localnodes_of_edges_in_cells(1, iiE, iK) = iii;
        }
      }
    }
    count++;
  }

  // arma::uvec badindices_nodes_of_edges = arma::find(nodes_of_edges=-1);
  // arma::uvec badindices_edges_of_cells = arma::find(edges_of_cells=-1);
  // arma::uvec badindices_localnodes = arma::find(localnodes_of_edges_in_cells=-1);
  // std::cerr << "badindices_nodes_of_edges " << badindices_nodes_of_edges.t();
  // std::cerr << "badindices_edges_of_cells " << badindices_edges_of_cells.t();
  // std::cerr << "badindices_localnodes " << badindices_localnodes.t();

  if(debug)
  {
    std::cerr << "nodes_of_cells\n";
    for(int iK=0; iK < ncells; iK++)
    {
      std::cerr << iK << " ==> " << alat::armaivec(nodes_of_cells.col(iK)).t();
    }
    std::cerr << "nodes_of_edges\n";
    for(int iE=0; iE < nedges; iE++)
    {
      std::cerr << iE << " ==> " << alat::armaivec(nodes_of_edges.col(iE)).t();
    }
    std::cerr << "edges_of_cells\n";
    for(int iK=0; iK < ncells; iK++)
    {
      std::cerr << iK << " ==> " << alat::armaivec(edges_of_cells.col(iK)).t();
      std::cerr <<  alat::armaimat(localnodes_of_edges_in_cells.slice(iK));
    }
  }
}

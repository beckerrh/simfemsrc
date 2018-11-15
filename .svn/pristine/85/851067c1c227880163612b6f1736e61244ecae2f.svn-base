#include  "Mesh/boundaryinformation.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Mesh/sidesconstructor.hpp"
#include  "Mesh/sidesandcells.hpp"
#include  <algorithm>
#include  <limits>

using namespace mesh;

/*---------------------------------------------------------*/
template<int NODESPERSIDE>
SidesConstructor<NODESPERSIDE>::~SidesConstructor(){}
template<int NODESPERSIDE>
SidesConstructor<NODESPERSIDE>::SidesConstructor(MeshUnitInterface* mesh) : _mesh(mesh) {}

/*---------------------------------------------------------*/
template<int NODESPERSIDE>
void SidesConstructor<NODESPERSIDE>::constructSidesFromCells(const MeshUnitInterface::BoundarySideToColorForSidesConstructor& bstcin, const alat::Map<int, alat::armaivec>* partition_to_cell, int color_default, bool debug)
{
  alat::Map<int, alat::IntSet> partition_to_cell_set;
  if(partition_to_cell)
  {
    for(alat::Map<int, alat::armaivec>::const_iterator p = partition_to_cell->begin(); p!=partition_to_cell->end();p++)
    {
      std::copy(p->second.begin(), p->second.end(), std::inserter(partition_to_cell_set[p->first], partition_to_cell_set[p->first].begin()));
    }
  }
  // std::cerr << "partition_to_cell_set = " << partition_to_cell_set << "\n";
  // std::cerr << "bstcin = " << bstcin << "\n";
  const alat::armaimat&  _nodes_of_cells = _mesh->getNodesAndNodesOfCells().getNodesOfCells();
  alat::armaimat&  _nodes_of_sides = _mesh->getSidesAndCells().getNodesOfSides();
  alat::armaimat&  _sides_of_cells = _mesh->getSidesAndCells().getSidesOfCells();
  alat::armaimat&  _cells_of_sides = _mesh->getSidesAndCells().getCellsOfSides();
  assert(NODESPERSIDE==_mesh->getNNodesPerSide());
  int SIDESPERCELL = _mesh->getNSidesPerCell();

  //! Assumption: no hanging nodes !!
  BoundarySideToColor bstc;
  for(MeshUnitInterface::BoundarySideToColorForSidesConstructor::const_iterator p=bstcin.begin();p!=bstcin.end();p++)
  {
    assert(p->second.n_rows==NODESPERSIDE);
    for(int i=0;i<p->second.n_cols;i++)
    {
      FixSide S(p->second.col(i));
      std::sort(S.begin(),S.end());
      bstc[S]=p->first;
    }
  }
  // std::cerr << "bstc = " << bstc << "\n";

  SideToInfo _found;

  // first round : how many interor sides ?
  // interior == in two cells
  // boundary == only in one cell
  int ninteriorsides = 0;
  for(int i = 0; i < _nodes_of_cells.n_cols; i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      FixSide s(_mesh->_getSideOfCell(i, ii));
      std::sort( s.begin(), s.end() );
      typename SideToInfo::iterator p = _found.find(s);
      if( p == _found.end() )
      {
        alat::FixArray<2, int> index;
        _found.insert( std::make_pair(s, index) );
      }
      else
      {
        ninteriorsides++;
        _found.erase(p);
      }
    }
  }
  int nbsides = _found.size();
  int nsides = ninteriorsides+nbsides;
  // std::cerr << "nbsides="<<nbsides<<" ninteriorsides="<<ninteriorsides<<"\n";
  _nodes_of_sides.set_size(NODESPERSIDE, nsides);

  // second round : insert data
  _found.clear();
  _sides_of_cells.set_size(SIDESPERCELL, _mesh->getNCells() );
  int count = 0;
  for(int i = 0; i < _nodes_of_cells.n_cols; i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      FixSide s(_mesh->_getSideOfCell(i, ii));
      FixSide ssort = s;
      std::sort( ssort.begin(), ssort.end() );
      typename SideToInfo::iterator p = _found.find(ssort);
      if( p == _found.end() )
      {
        alat::FixArray<2, int> index;
        index[0] = i;
        index[1] = ii;
        _found.insert( std::make_pair(ssort, index) );
      }
      else
      {
        // we put in the NON-SORTED !
        for(int titi=0;titi<NODESPERSIDE;titi++) _nodes_of_sides(titi,count) = s[titi];
        int k = p->second[0];
        int kk = p->second[1];
        _sides_of_cells(kk,k) = count;
        _sides_of_cells(ii,i) = count;
        _found.erase(p);
        count++;
      }
    }
  }
  // now _found contains all boundarysides (sides found once)
  // we add all boundary-sides, which are not yet defined, by giving it the color color_default !!
  int nadditionalboundarysides = 0;
  for(int i = 0; i < _nodes_of_cells.n_cols; i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      FixSide ssort(_mesh->_getSideOfCell(i, ii));
      std::sort( ssort.begin(), ssort.end() );
      typename SideToInfo::iterator p = _found.find(ssort);
      if( p == _found.end() ) {continue;}
      // std::cerr << "i ii " << i << " " << ii << " p->second " << p->second[0] << " " << p->second[1] << "\n";
      assert(i == p->second[0]);
      assert(ii == p->second[1]);
      typename BoundarySideToColor::const_iterator pb = bstc.find(ssort);
      if( pb != bstc.end() ) {continue;}
      bool foundinparts = false;
      for(alat::Map<int, alat::IntSet>::const_iterator q=partition_to_cell_set.begin();q!=partition_to_cell_set.end();q++)
      {
        if(q->second.find(i)!=q->second.end())
        {
          bstc[ssort] =  std::numeric_limits<int>::lowest();
          foundinparts=true;
          break;
        }
      }
      if(not foundinparts)
      {
        if(debug) std::cerr<<"*** not found "<<ssort<< " for cell " << i << "\n";
        // std::cerr << "bstc = " << bstc << "\n";
        // assert(0);
        bstc[ssort] = color_default;
        nadditionalboundarysides++;
      }
    }
  }
  if(nadditionalboundarysides)
  {
    std::cout<<"constructSidesFromCells(): # boundary sides without color (given color " <<  color_default << " ): "<<nadditionalboundarysides<<std::endl;
  }
  // std::cerr << "bstc = " << bstc << "\n";

  alat::IntMap size_of_color;
  // for(typename BoundarySideToColor::const_iterator p = bstc.begin(); p != bstc.end(); p++)
  // {
  //   int col = p->second;
  //   if( size_of_color.find(col) == size_of_color.end() )
  //   {
  //     size_of_color[col] = 0;
  //   }
  //   size_of_color[col]++;
  // }
  for(int i = 0; i < _nodes_of_cells.n_cols; i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      FixSide s(_mesh->_getSideOfCell(i, ii));
      FixSide ssort = s;
      std::sort( ssort.begin(), ssort.end() );
      typename SideToInfo::iterator p = _found.find(ssort);
      if( p == _found.end() ) {continue;}
      typename BoundarySideToColor::const_iterator pb = bstc.find(p->first);
      assert( pb != bstc.end() );
      int color = pb->second;
      if(color==std::numeric_limits<int>::lowest())
      {
        for(alat::Map<int, alat::IntSet>::const_iterator q=partition_to_cell_set.begin();q!=partition_to_cell_set.end();q++)
        {
          if(q->second.find(i)!=q->second.end())
          {
            if( size_of_color.find(q->first) == size_of_color.end() )
            {
              size_of_color[q->first] = 0;
            }
            size_of_color[q->first]++;
          }
        }
      }
      else
      {
        if( size_of_color.find(color) == size_of_color.end() )
        {
          size_of_color[color] = 0;
        }
        size_of_color[color]++;
      }
    }
  }
  // std::cerr << "size_of_color = " << size_of_color << "\n";
  alat::IntMap col2size_bsides;
  for(std::map<int, int>::const_iterator p = size_of_color.begin();p!=size_of_color.end();p++)
  {
    if(_mesh->getBoundaryInformation(p->first).size()==0)
    {
      _mesh->getBoundaryInformation(p->first).set_size(p->second);
    }
    assert(p->second==_mesh->getBoundaryInformation(p->first).size());
    col2size_bsides[p->first] = 0;
  }
  // insert _nodes_of_sides and _sides_of_cells four bdry sides
  for(int i = 0; i < _nodes_of_cells.n_cols; i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      FixSide s(_mesh->_getSideOfCell(i, ii));
      FixSide ssort = s;
      std::sort( ssort.begin(), ssort.end() );
      typename SideToInfo::iterator p = _found.find(ssort);
      if( p == _found.end() ) {continue;}
      typename BoundarySideToColor::const_iterator pb = bstc.find(p->first);
      assert( pb != bstc.end() );
      for(int titi=0;titi<NODESPERSIDE;titi++) _nodes_of_sides(titi,count) = s[titi];
      assert(i == p->second[0]);
      assert(ii == p->second[1]);
      _sides_of_cells(ii,i) = count;
      count++;
    }
  }
  for(int i = 0; i < _nodes_of_cells.n_cols; i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      FixSide s(_mesh->_getSideOfCell(i, ii));
      FixSide ssort = s;
      std::sort( ssort.begin(), ssort.end() );
      typename SideToInfo::iterator p = _found.find(ssort);
      if( p == _found.end() ) {continue;}
      typename BoundarySideToColor::const_iterator pb = bstc.find(p->first);
      assert( pb != bstc.end() );
      int color = pb->second;
      if(color==std::numeric_limits<int>::lowest())
      {
        for(alat::Map<int, alat::IntSet>::const_iterator q=partition_to_cell_set.begin();q!=partition_to_cell_set.end();q++)
        {
          if(q->second.find(i)!=q->second.end())
          {
            alat::armaimat& cells_on_bdry = _mesh->getBoundaryInformation(q->first).getCellsOnBdryOfPlain();
            assert(i == p->second[0]);
            assert( ii == p->second[1]);
            cells_on_bdry(0, col2size_bsides[q->first]) = i;
            cells_on_bdry(1, col2size_bsides[q->first]) = _sides_of_cells(ii,i);
            cells_on_bdry(2, col2size_bsides[q->first]) = ii;
            col2size_bsides[q->first]++;
          }
        }
      }
      else
      {
        alat::armaimat& cells_on_bdry = _mesh->getBoundaryInformation(color).getCellsOnBdryOfPlain();
        // std::cerr << "boundary side found " << s << "\n";
        assert(i == p->second[0]);
        assert( ii == p->second[1]);
        cells_on_bdry(0, col2size_bsides[color]) = i;
        cells_on_bdry(1, col2size_bsides[color]) = _sides_of_cells(ii,i);
        cells_on_bdry(2, col2size_bsides[color]) = ii;
        col2size_bsides[color]++;
      }
    }
  }
  _cells_of_sides.set_size(2,nsides);
  for(int is = 0; is < nsides; is++)
  {
    _cells_of_sides(0,is) = -1;
    _cells_of_sides(1,is) = -1;
  }
  for(int i = 0; i < _nodes_of_cells.n_cols; i++)
  {
    for(int ii = 0; ii < SIDESPERCELL; ii++)
    {
      int is = _sides_of_cells(ii,i);
      if(_cells_of_sides(0,is) == -1)
      {
        _cells_of_sides(0,is) = i;
      }
      else
      {
        int help = _cells_of_sides(0,is);
        _cells_of_sides(0,is) = i;
        _cells_of_sides(1,is) = help;
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
// line mesh
#define NODESPERSIDE 1
template class mesh::SidesConstructor<NODESPERSIDE>;
#undef NODESPERSIDE

// triangle mesh
#define NODESPERSIDE 2
template class mesh::SidesConstructor<NODESPERSIDE>;
#undef NODESPERSIDE

// tetrahedral mesh
#define NODESPERSIDE 3
template class mesh::SidesConstructor<NODESPERSIDE>;
#undef NODESPERSIDE

// hexahedral mesh
#define NODESPERSIDE 4
template class mesh::SidesConstructor<NODESPERSIDE>;
#undef NODESPERSIDE

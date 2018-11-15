#include  "FadalightAdaptiveMesh/multimeshcreator.hpp"
#include  "FadalightMesh/refineinfo.hpp"
#include  <sstream>

using namespace FadalightAdaptiveMesh;
using namespace std;

/*--------------------------------------------------------------------------*/
MultiMeshCreator::~MultiMeshCreator() {}
MultiMeshCreator::MultiMeshCreator(Coarsener2d* coarsener): _coarsener(coarsener){}
MultiMeshCreator::MultiMeshCreator(const MultiMeshCreator& multimeshcreator) {assert(0);}
MultiMeshCreator& MultiMeshCreator::operator= (const MultiMeshCreator & multimeshcreator) {assert(0);}
std::string MultiMeshCreator::getClassName() const
{
    return "MultiMeshCreator";
}

AdaptiveMeshInterface* MultiMeshCreator::_getMesh(){return _coarsener->getMesh();}
tree<FaceInterface*> & MultiMeshCreator::_getFaces(){return _coarsener->getFaces();}

/*--------------------------------------------------------------------------*/
void MultiMeshCreator::createMultiMesh(std::string dirname, int nlevels, int ncells, arma::file_type datatype)
{
  string multimeshname = dirname+".fadalightmesh/MultiMesh";
  string cmd = "rm -rf "+multimeshname;
  int error = system( cmd.c_str() );
  assert(!error);
  cmd = "mkdir "+multimeshname;
  error = system( cmd.c_str() );
  assert(!error);

  int count = 1;
  // bool coarsen = 1;
  _getMesh()->createGeometryObject("RefineInfo");
  FadalightMesh::RefineInfo* _refineinfo= dynamic_cast<FadalightMesh::RefineInfo*>(_getMesh()->getGeometryObject("RefineInfo"));
  // if(_getFaces().max_depth() == 0)
  // {
  //   coarsen = 0;
  // }
  // while(coarsen)
  while(( _getFaces().max_depth() > 0 ) && (count < nlevels)  && (_getFaces().size() > ncells))
  {
    _coarsener->globalCoarsen();
    _constructRefineInfo();
    std::stringstream ss0;
    ss0<<count;
    std::string outfilename = multimeshname+"/mesh."+ss0.str();
    // coarsen = ( _getFaces().max_depth() > 0 ) && (count < nlevels)  && (_getFaces().size() > ncells);
    _getMesh()->writeFadalightMesh(outfilename, datatype);
    // if(count == 1)
    // {
    //   outfilename = dirname+".fadalightmesh/RefineInfoML";
    //   _refineinfo->write(outfilename, datatype);
    // }
    // else
    // {
    //   std::stringstream ss1;
    //   ss1<<count-1;
    //   outfilename = multimeshname+"/mesh."+ss1.str()+".fadalightmesh/RefineInfoML";
    //   _refineinfo->write(outfilename, datatype);
    // }
    count++;
  }

  string filename = multimeshname+"/n";
  ofstream file( filename.c_str() );
  assert( file.is_open() );
  file<<count<<"\n";
  file.close();
}

/*--------------------------------------------------------------------------*/

void MultiMeshCreator::_constructRefineInfo()
{
  int n = _coarsener->getFaceCoarse().size();
  alat::Vector<alat::armaivec> SPC(n), SPN(n), SPE(n);
  for(face_leafpointer f = _getFaces().begin_leaf(); f != _getFaces().end_leaf(); f++)
  {
    int i =  _coarsener->getFaceId2Id()[( *f )->id()];

    SPC[i].set_size(_coarsener->getFaceCoarse()[f].size());
    SPC[i] = _coarsener->getFaceCoarse()[f];

    SPN[i].set_size(_coarsener->getNodeCoarse()[f].size());
    SPN[i] = _coarsener->getNodeCoarse()[f];

    SPE[i].set_size(_coarsener->getEdgeCoarse()[f].size());
    SPE[i] = _coarsener->getEdgeCoarse()[f];
  }
  FadalightMesh::RefineInfo* _refineinfo= dynamic_cast<FadalightMesh::RefineInfo*>(_getMesh()->getGeometryObject("RefineInfo"));
  alat::SparsityPattern& coarsenodeids = _refineinfo->getCoarseNodes();
  alat::SparsityPattern& coarsecellids = _refineinfo->getCoarseCells();
  alat::SparsityPattern& coarsesideids = _refineinfo->getCoarseSides();

  coarsenodeids.set_size(SPN);
  coarsecellids.set_size(SPC);
  coarsesideids.set_size(SPE);

  alat::armaivec& nodeids = _refineinfo->getNodeIds();
  alat::armaivec& sideids = _refineinfo->getSideIds();
  alat::armaivec& cellids = _refineinfo->getCellIds();
  const alat::IntMap& nodeid2id = _coarsener->_getNodeId2Id();
  const alat::IntMap& edgeid2id = _coarsener->getEdgeId2Id();
  const alat::IntMap& faceid2id = _coarsener->getFaceId2Id();

  // std::cerr << "nodeid2id="<<nodeid2id<<"\n";
  // std::cerr << "edgeid2id="<<edgeid2id<<"\n";
  // std::cerr << "faceid2id="<<faceid2id<<"\n";

  nodeids.resize(nodeid2id.size());
  for(alat::IntMap::const_iterator p= nodeid2id.begin(); p!= nodeid2id.end(); p++)
  {
    nodeids[p->second] = p->first;
  }
  sideids.resize(edgeid2id.size());
  for(alat::IntMap::const_iterator p= edgeid2id.begin(); p!= edgeid2id.end(); p++)
  {
    sideids[p->second] = p->first;
  }
  cellids.resize(faceid2id.size());
  for(alat::IntMap::const_iterator p= faceid2id.begin(); p!= faceid2id.end(); p++)
  {
    cellids[p->second] = p->first;
  }
}

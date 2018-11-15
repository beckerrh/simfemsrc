#include  "FadalightAdaptiveMesh/multimeshcreator3d.hpp"
#include  "FadalightAdaptiveMesh/volumeinterface.hpp"
#include  "FadalightMesh/refineinfo.hpp"
#include  <sstream>

using namespace FadalightAdaptiveMesh;
using namespace std;

/*--------------------------------------------------------------------------*/
void MultiMeshCreator3d::createMultiMesh(std::string dirname, int nlevels, int ncells, arma::file_type datatype)
{
  string multimeshname = dirname+".fadalightmesh/MultiMesh";
  string cmd = "rm -rf "+multimeshname;
  int error = system( cmd.c_str() );
  assert(!error);
  cmd = "mkdir "+multimeshname;
  error = system( cmd.c_str() );
  assert(!error);

  int count = 1;
  bool coarsen = 1;
  _getMesh()->createGeometryObject("RefineInfo");
  FadalightMesh::RefineInfo* _refineinfo= dynamic_cast<FadalightMesh::RefineInfo*>(_getMesh()->getGeometryObject("RefineInfo"));
  if(_getVolumes().max_depth() == 0)
  {
    coarsen = 0;
  }
  while(coarsen)
  {
    _coarsener->globalCoarsen();
    _constructRefineInfo();
    std::stringstream ss0;
    ss0<<count;
    std::string outfilename = multimeshname+"/mesh."+ss0.str();
    coarsen = ( _getVolumes().max_depth() > 0 );
    _getMesh()->writeFadalightMesh(outfilename, datatype);
    if(count == 1)
    {
      outfilename = dirname+".fadalightmesh/RefineInfo";
      _refineinfo->save(outfilename, datatype);
    }
    else
    {
      std::stringstream ss1;
      ss1<<count-1;
      outfilename = multimeshname+"/mesh."+ss1.str()+".fadalightmesh/RefineInfo";
      _refineinfo->save(outfilename, datatype);
    }
    count++;
  }

  string filename = multimeshname+"/n";
  ofstream file( filename.c_str() );
  assert( file.is_open() );
  file<<count<<"\n";
  file.close();
}

/*--------------------------------------------------------------------------*/

void MultiMeshCreator3d::_constructRefineInfo()
{
  int n = _coarsener->getVolumeCoarse().size();
  alat::Vector<alat::armaivec> SPC(n), SPF(n), SPE(n), SPN(n) ;
  for(volume_leafpointer v = _getVolumes().begin_leaf(); v != _getVolumes().end_leaf(); v++)
  {
    int i =  _coarsener->getVolumeId2Id()[( *v )->id()];

    SPC[i].set_size(_coarsener->getVolumeCoarse()[v].size());
    SPC[i] = _coarsener->getVolumeCoarse()[v];

    SPF[i].set_size(_coarsener->getFaceCoarse()[v].size());
    SPF[i] = _coarsener->getFaceCoarse()[v];

    SPE[i].set_size(_coarsener->getEdgeCoarse()[v].size());
    SPE[i] = _coarsener->getEdgeCoarse()[v];

    SPN[i].set_size(_coarsener->getNodeCoarse()[v].size());
    SPN[i] = _coarsener->getNodeCoarse()[v];


  }
  FadalightMesh::RefineInfo* _refineinfo= dynamic_cast<FadalightMesh::RefineInfo*>(_getMesh()->getGeometryObject("RefineInfo"));

  alat::SparsityPattern& coarsenodeids = _refineinfo->getCoarseNodes();
  alat::SparsityPattern& coarseedgeids = _refineinfo->getCoarseEdges();
  alat::SparsityPattern& coarsesideids = _refineinfo->getCoarseSides();
  alat::SparsityPattern& coarsecellids = _refineinfo->getCoarseCells();
  coarsenodeids.set_size(SPN);
  coarseedgeids.set_size(SPE);
  coarsesideids.set_size(SPF);
  coarsecellids.set_size(SPC);
}

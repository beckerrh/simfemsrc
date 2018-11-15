#include  "Alat/vectoronevariable.hpp"
#include  "Mesh/meshvisitorinterface.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Mesh/measureofcell.hpp"
#include  "Mesh/nodesandnodesofcells.hpp"
#include  "Mesh/enums.hpp"
#include  "Mesh/sidesandcells.hpp"
#include  "Mesh/normals.hpp"
#include  <cassert>
#include  <mpi.h>

using namespace mesh;

/*--------------------------------------------------------------------------*/
CutInterface::~CutInterface(){}
CutInterface::CutInterface() : GeometryObject()
{
  // std::cerr << "CutInterface::CutInterface()\n";
}
CutInterface::CutInterface( const CutInterface& geometryobject) : GeometryObject(geometryobject)
{
  std::cerr << "CutInterface::CutInterface(geometryobject)\n";
  _cutcells = geometryobject._cutcells;
  _cutedges = geometryobject._cutedges;
  _celliscut = geometryobject._celliscut;
  _edgeiscut = geometryobject._edgeiscut;
  _cutnodes = geometryobject._cutnodes;
  _cutcoeff = geometryobject._cutcoeff;
  _nodesofcutcellsisin = geometryobject._nodesofcutcellsisin;
  _nodesofcutcell = geometryobject._nodesofcutcell;
  _measuresofcutcells = geometryobject._measuresofcutcells;
  _normalsofcutcells = geometryobject._normalsofcutcells;
}
CutInterface& CutInterface::operator=( const CutInterface& geometryobject)
{
  std::cerr << "CutInterface::operator=(geometryobject)\n";
  GeometryObject::operator=(geometryobject);
  assert(0);
  return *this;
}
std::string CutInterface::getClassName() const
{
  return "CutInterface";
}
std::unique_ptr<GeometryObject> CutInterface::clone() const
{
  std::cerr << "CutInterface::clone()\n";
  assert(0);
  return std::unique_ptr<mesh::GeometryObject>(new CutInterface(*this));
}

/*--------------------------------------------------------------------------*/
const alat::armaivec& CutInterface::getCutCells() const {return _cutcells;}
const alat::armaivec& CutInterface::getCutEdges() const {return _cutedges;}
const alat::armaivec& CutInterface::getCellIsCut() const {return _celliscut;}
const alat::armaivec& CutInterface::getEdgeIsCut() const {return _edgeiscut;}
const alat::armaivec& CutInterface::getNodeIsCut() const {return _nodeiscut;}
const alat::armaivec& CutInterface::getCutNodeIsIn() const {return _cutnodesisin;}
const alat::armaivec& CutInterface::getCutNodes() const {return _cutnodes;}
const arma::vec& CutInterface::getCutCoeff() const {return _cutcoeff;}
const alat::armaimat& CutInterface::getNodesOfCutCellsIsIn() const {return _nodesofcutcellsisin;}
const alat::armaimat& CutInterface::getNodesOfCutCells() const {return _nodesofcutcell;}
const arma::mat& CutInterface::getMeasureOfCutCells() const {return _measuresofcutcells;}
const arma::mat& CutInterface::getNormalsOfCutCells() const {return _normalsofcutcells;}
const arma::mat& CutInterface::getCInOfCutCells() const {return _cinofcutcells;}
const arma::mat& CutInterface::getCExOfCutCells() const {return _cexofcutcells;}
const arma::mat& CutInterface::getCofInOfCutCells() const {return _cofinofcutcells;}
const arma::mat& CutInterface::getCofExOfCutCells() const {return _cofexofcutcells;}
/*--------------------------------------------------------------------------*/
alat::armaivec CutInterface::getSizes() const
{
  alat::armaivec sizes(1);
  sizes[0] = _cutcells.n_elem;
  return sizes;
}
void CutInterface::setSizes(alat::armaivec::const_iterator sizes)
{
  _cutcells.set_size(sizes[0]);
}
void CutInterface::send(int neighbor, int tag) const
{
  MPI_Request request;
  MPI_Isend(_cutcells.begin(), _cutcells.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
}
void CutInterface::recv(int neighbor, int tag)
{
  MPI_Status status;
  MPI_Request request;
  MPI_Irecv(_cutcells.begin(), _cutcells.size(), MPI_INT, neighbor, tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  // std::cerr << "CutInterface::recv() " << _normals.t() << "\n";
}

/*--------------------------------------------------------------------------*/
void CutInterface::loadH5(const arma::hdf5_name& spec)
{
  _cutcells.load(arma::hdf5_name(spec.filename, spec.dsname+"/cutcells", spec.opts));
}
void CutInterface::saveH5(const arma::hdf5_name& spec) const
{
  _cutcells.save(arma::hdf5_name(spec.filename, spec.dsname+"/cutcells", spec.opts));
}

/*--------------------------------------------------------------------------*/
double CutInterface::_surface(arma::subview_col<double> u, arma::subview_col<double> v, arma::subview_col<double> w)const
{
  arma::vec p(3), q(3), normal(3);
  p[0] = v[0]-u[0];
  p[1] = v[1]-u[1];
  p[2] = v[2]-u[2];
  q[0] = w[0]-u[0];
  q[1] = w[1]-u[1];
  q[2] = w[2]-u[2];
  normal[0] = p[1]*q[2]-p[2]*q[1];
  normal[1] = p[2]*q[0]-p[0]*q[2];
  normal[2] = p[0]*q[1]-p[1]*q[0];
  return arma::norm(normal)*0.5;
}
/*--------------------------------------------------------------------------*/
double CutInterface::_surface(const arma::mat& nodes)const
{
  assert(nodes.n_cols==4);
  double vol = 0.0;
  vol += _surface(nodes.col(1),nodes.col(2),nodes.col(3));
  vol += _surface(nodes.col(0),nodes.col(2),nodes.col(3));
  vol += _surface(nodes.col(0),nodes.col(1),nodes.col(3));
  vol += _surface(nodes.col(0),nodes.col(1),nodes.col(2));
  return 0.5*vol;
}
/*--------------------------------------------------------------------------*/
void CutInterface::_computeCutTet(double& volin, double& volex, arma::subview_col<double> xcin, arma::subview_col<double> xcex, const arma::vec& normal, int iK, const mesh::MeshUnitInterface* mesh, const arma::mat& innodes, const arma::mat& exnodes, const arma::mat& addnodes, const arma::ivec& innodesind, const arma::ivec& exnodesind, const arma::ivec& addnodesedges)const
{
  const SidesAndCells& sidesandcells=mesh->getSidesAndCells();
  const alat::armaimat& nodesofsides = sidesandcells.getNodesOfSides();
  const alat::armaimat& sidesofcells = sidesandcells.getSidesOfCells();
  const EdgesAndCells& edgesandcells = mesh->getEdgesAndCells();
  const alat::armaimat& nodesofedges = edgesandcells.getNodesOfEdges();
  const NodesAndNodesOfCells& nodesandnodesofcells = mesh->getNodesAndNodesOfCells();
  const alat::armaimat& nodesofcells = nodesandnodesofcells.getNodesOfCells();
  const arma::mat& nodes = nodesandnodesofcells.getNodes();
  assert(mesh->geometryObjectExists(meshEnums::Normals));
  std::shared_ptr<const mesh::Normals> normalsobj = std::dynamic_pointer_cast<const mesh::Normals>(mesh->getGeometryObject(meshEnums::Normals));
  assert(normalsobj);
  const arma::mat& normals = normalsobj->getNormals();
  const arma::fmat& sigmas = normalsobj->getSigma();

  assert(sidesofcells.n_rows==4);
  assert(nodesofsides.n_rows==3);
  assert(sigmas.n_rows==4);

  alat::IntMap innodeindtoind, exnodeindtoind, raplacein, replaceex;
  assert(nodesofcells.n_rows==4);
  for(int ii=0;ii<innodesind.size();ii++)
  {
    innodeindtoind[innodesind[ii]]=ii;
  }
  for(int ii=0;ii<exnodesind.size();ii++)
  {
    exnodeindtoind[exnodesind[ii]]=ii;
  }
  for(int ii=0;ii<addnodesedges.size();ii++)
  {
    int iE = addnodesedges[ii];
    int iN0 = nodesofedges(0,iE);
    int iN1 = nodesofedges(1,iE);
    if(innodeindtoind.find(iN0)==innodeindtoind.end())
    {
      raplacein[iN0] = ii;
      replaceex[iN1] = ii;
    }
    else
    {
      raplacein[iN1] = ii;
      replaceex[iN0] = ii;
    }
  }

  double vol=0.0;
  volin = 0.0;
  volex = 0.0;
  xcin.zeros();
  xcex.zeros();
  arma::vec xcs(3), xcsin(3), xcsex(3), xK(3), xK2(3);
  xK.zeros();
  for(int ii=0;ii<4;ii++)
  {
    int iN = nodesofcells(ii,iK);
    xK +=nodes.col(iN)/4.0;
  }
  xK2.zeros();
  for(int iis=0;iis<4;iis++)
  {
    int iS = sidesofcells(iis, iK);
    arma::subview_col<double> normal = normals.col(iS);
    double sigma = sigmas(iis,iK);
    double dnorm = arma::norm(normal);
    int iN = nodesofsides(0,iS);
    double scale = sigma*arma::dot(nodes.col(iN), normal)/3.0/dnorm;
    vol += dnorm*scale;

    arma::mat snodes(3,3);
    double surfin, surfex;
    int count=0;
    xcs.zeros();
    for(int ii=0;ii<3;ii++)
    {
      int iN = nodesofsides(ii,iS);
      xcs += nodes.col(iN)/3.0;
      if(innodeindtoind.find(iN)!=innodeindtoind.end()) count++;
    }
    xK2 += dnorm*scale*xcs*3.0/4.0;
    if(count==1)
    {
      for(int ii=0;ii<3;ii++)
      {
        int iN = nodesofsides(ii,iS);
        if(innodeindtoind.find(iN)!=innodeindtoind.end()) snodes.col(ii)=nodes.col(iN);
        else snodes.col(ii)=addnodes.col(raplacein[iN]);
      }
      surfin = _surface(snodes.col(0), snodes.col(1), snodes.col(2));
      assert(surfin<dnorm);
      surfex = dnorm - surfin;
      xcsin = (snodes.col(0) + snodes.col(1) + snodes.col(2))/3.0;
      xcsex = (dnorm/surfex)*xcs -(surfin/surfex)*xcsin;
    }
    else
    {
      for(int ii=0;ii<3;ii++)
      {
        int iN = nodesofsides(ii,iS);
        if(exnodeindtoind.find(iN)!=exnodeindtoind.end()) snodes.col(ii)=nodes.col(iN);
        else snodes.col(ii)=addnodes.col(replaceex[iN]);
      }
      surfex = _surface(snodes.col(0), snodes.col(1), snodes.col(2));
      assert(surfex<dnorm);
      surfin = dnorm - surfex;
      xcsex = (snodes.col(0) + snodes.col(1) + snodes.col(2))/3.0;
      xcsin = (dnorm/surfin)*xcs -(surfex/surfin)*xcsex;
    }
    volin += scale*surfin;
    volex += scale*surfex;
    xcin += scale*surfin*xcsin*3.0/4.0;
    xcex += scale*surfex*xcsex*3.0/4.0;
  }
  xK2/= vol;
  xcin/= volin;
  xcex/= volex;
  if(arma::norm(xK-xK2)>1e-12)
  {
    std::cerr << "xK = " << xK.t() << " xK2 = " << xK2.t();
    assert(0);
  }

  double d = arma::dot(normal,addnodes.col(0))/3.0;
  volin += d;
  volex -= d;
  xcs = (addnodes.col(0) + addnodes.col(1) + addnodes.col(2) + addnodes.col(3))/4.0;
  xcin += d*xcs*3.0/4.0;
  xcex -= d*xcs*3.0/4.0;

  if(fabs(vol-volin-volex)>1e-12 or volin <=0.0 or volex <=0.0)
  {
    std::cerr << "vol = " << vol << " volin = " << volin << " volex = " << volex << "\n";
    std::cerr << "d = " << d << "\n";
    std::cerr << "vol-volin-volex = " << vol-volin-volex << "\n";
    // assert(0);
  }
  if( arma::norm(xK - (volin/vol)*xcin - (volex/vol)*xcex ) > 1e-12)
  {
    std::cerr << "xK="<<xK.t();
    arma::vec xt = (volin/vol)*xcin + (volex/vol)*xcex;
    std::cerr << "xt="<<xt.t();
    std::cerr << "xcin="<<xcin.t();
    std::cerr << "xcex="<<xcex.t();
    // assert(0);
  }
}
/*--------------------------------------------------------------------------*/
void CutInterface::_computeMeasuresOfCutCell(int iK, const mesh::MeshUnitInterface* mesh, arma::subview_col<double> measuresofcutcell, arma::subview_col<double> normalsofcutcells, arma::subview_col<double> cinofcutcells, arma::subview_col<double> cexofcutcells, const arma::mat&  innodes, const arma::mat& exnodes, const arma::mat& addnodes, double moc, const arma::ivec& innodesind, const arma::ivec& exnodesind, const arma::ivec& addnodesedges)
{
  if(mesh->getDimension()==2)
  {
    assert(innodes.n_cols+exnodes.n_cols==3);
    assert(addnodes.n_cols==2);
    arma::vec normal(3), test(3);
    normal[0] = addnodes(1,1)-addnodes(1,0);
    normal[1] = addnodes(0,0)-addnodes(0,1);
    normal[2] = 0.0;
    test[0] = 0.5*(addnodes(0,0)+addnodes(0,1));
    test[1] = 0.5*(addnodes(1,0)+addnodes(1,1));
    test[2] = 0.0;
    arma::vec xK(3, arma::fill::zeros);
    for(int ii=0;ii<innodes.n_cols;ii++)
    {
      xK += innodes.col(ii)/3.0;
    }
    for(int ii=0;ii<exnodes.n_cols;ii++)
    {
      xK += exnodes.col(ii)/3.0;
    }
    arma::mat nodes(2,2);
    if(innodes.n_cols==1)
    {
      for(int ii=0;ii<2;ii++)
      {
        nodes(0,ii) = addnodes(0,ii)-innodes(0,0);
        nodes(1,ii) = addnodes(1,ii)-innodes(1,0);
      }
      measuresofcutcell(0) = fabs(arma::det(nodes))/2.0;
      measuresofcutcell(1) = moc - measuresofcutcell(0);
      test[0] -= innodes(0,0);
      test[1] -= innodes(1,0);
      if(arma::dot(test, normal)<0.0)
      {
        normal *= -1.0;
      }
      cinofcutcells = (1.0/3.0)*(innodes.col(0) + addnodes.col(0) + addnodes.col(1));
      cexofcutcells = (moc/measuresofcutcell(1))*xK - (measuresofcutcell(0)/measuresofcutcell(1))*cinofcutcells;
    }
    else
    {
      assert(exnodes.n_cols==1);
      for(int ii=0;ii<2;ii++)
      {
        nodes(0,ii) = addnodes(0,ii)-exnodes(0,0);
        nodes(1,ii) = addnodes(1,ii)-exnodes(1,0);
      }
      measuresofcutcell(1) = fabs(arma::det(nodes))/2.0;
      measuresofcutcell(0) = moc - measuresofcutcell(1);
      test[0] -= exnodes(0,0);
      test[1] -= exnodes(1,0);
      if(arma::dot(test, normal)>0.0)
      {
        normal *= -1.0;
      }
      cexofcutcells = (1.0/3.0)*(exnodes.col(0) + addnodes.col(0) + addnodes.col(1));
      cinofcutcells = (moc/measuresofcutcell(0))*xK - (measuresofcutcell(1)/measuresofcutcell(0))*cexofcutcells;
    }
    normalsofcutcells = normal;
  }
  else if(mesh->getDimension()==3)
  {
    assert(addnodes.n_cols>=3);
    // std::cerr << "innodes = " << innodes;
    // std::cerr << "exnodes = " << exnodes;
    // std::cerr << "addnodes = " << addnodes;
    assert(innodes.n_cols+exnodes.n_cols==4);
    arma::vec xK(3, arma::fill::zeros);
    for(int ii=0;ii<innodes.n_cols;ii++)
    {
      xK += innodes.col(ii)/4.0;
    }
    for(int ii=0;ii<exnodes.n_cols;ii++)
    {
      xK += exnodes.col(ii)/4.0;
    }
    arma::vec p(3), q(3), normal(3), test(3);
    p[0] = addnodes(0,1)-addnodes(0,0);
    p[1] = addnodes(1,1)-addnodes(1,0);
    p[2] = addnodes(2,1)-addnodes(2,0);
    q[0] = addnodes(0,2)-addnodes(0,0);
    q[1] = addnodes(1,2)-addnodes(1,0);
    q[2] = addnodes(2,2)-addnodes(2,0);
    normal[0] = p[1]*q[2]-p[2]*q[1];
    normal[1] = p[2]*q[0]-p[0]*q[2];
    normal[2] = p[0]*q[1]-p[1]*q[0];
    normal *= 0.5;
    test.zeros();
    for(int ii=0;ii<addnodes.n_cols;ii++)
    {
      test += addnodes.col(ii)/(double)addnodes.n_cols;
    }
    if(innodes.n_cols==1)
    {
      arma::mat nodes(3,3);
      test -= innodes.col(0);
      for(int ii=0;ii<3;ii++)
      {
        nodes.col(ii) = addnodes.col(ii)-innodes.col(0);
      }
      measuresofcutcell(0) = fabs(arma::det(nodes))/6.0;
      measuresofcutcell(1) = moc - measuresofcutcell(0);
      if(arma::dot(test, normal)<0.0)
      {
        normal *= -1.0;
      }
      cinofcutcells = (1.0/4.0)*(innodes.col(0) + addnodes.col(0) + addnodes.col(1) + addnodes.col(2));
      cexofcutcells = (moc/measuresofcutcell(1))*xK - (measuresofcutcell(0)/measuresofcutcell(1))*cinofcutcells;
    }
    else if(exnodes.n_cols==1)
    {
      arma::mat nodes(3,3);
      test -= exnodes.col(0);
      for(int ii=0;ii<3;ii++)
      {
        nodes.col(ii) = addnodes.col(ii)-exnodes.col(0);
      }
      measuresofcutcell(1) = fabs(arma::det(nodes))/6.0;
      measuresofcutcell(0) = moc - measuresofcutcell(1);
      if(arma::dot(test, normal)>0.0)
      {
        normal *= -1.0;
      }
      cexofcutcells = (1.0/4.0)*(exnodes.col(0) + addnodes.col(0) + addnodes.col(1) + addnodes.col(2));
      cinofcutcells = (moc/measuresofcutcell(0))*xK - (measuresofcutcell(1)/measuresofcutcell(0))*cexofcutcells;
    }
    else
    {
      assert(innodes.n_cols==2);
      assert(exnodes.n_cols==2);
      assert(addnodes.n_cols==4);
      // correct normal length
      double surf = _surface(addnodes);
      double norm = arma::norm(normal);
      assert(surf>norm);
      normal *= surf/norm;
      test -=  0.5*(innodes.col(0) + innodes.col(1));
      if(arma::dot(test, normal)<0.0)
      {
        normal *= -1.0;
      }
      _computeCutTet(measuresofcutcell(0), measuresofcutcell(1), cinofcutcells, cexofcutcells, normal, iK, mesh, innodes, exnodes, addnodes, innodesind, exnodesind, addnodesedges);
    }
    normalsofcutcells = normal;
  }
  else
  {
    assert(0);
  }
}

/*--------------------------------------------------------------------------*/
void CutInterface::construct(const mesh::MeshUnitInterface* mesh, const GeometryConstructorInterface* geometryconstructor)
{
  // arma::mat testnodes(3,4);
  // testnodes(0,0) = 0.0; testnodes(1,0) = 0.0; testnodes(2,0) = 0.0;
  // testnodes(0,1) = 0.0; testnodes(1,1) = 1.0; testnodes(2,1) = 0.0;
  // testnodes(0,2) = 1.0; testnodes(1,2) = 1.0; testnodes(2,2) = 0.0;
  // testnodes(0,3) = 1.0; testnodes(1,3) = 0.0; testnodes(2,3) = 0.0;
  // std::cerr << "CutInterface::construct() _surfacetri="<<_surface(testnodes.col(1), testnodes.col(2),testnodes.col(3))<<"\n";
  // std::cerr << "CutInterface::construct() _surfacetri="<<_surface(testnodes.col(0), testnodes.col(2),testnodes.col(3))<<"\n";
  // std::cerr << "CutInterface::construct() _surfacetri="<<_surface(testnodes.col(0), testnodes.col(1),testnodes.col(3))<<"\n";
  // std::cerr << "CutInterface::construct() _surfacetri="<<_surface(testnodes.col(0), testnodes.col(1),testnodes.col(2))<<"\n";
  // std::cerr << "CutInterface::construct() _surfacequad="<<_surface(testnodes)<<"\n";
  // assert(0);



  // std::cerr << "CutInterface::construct() geometryconstructor="<< geometryconstructor->getClassName()<<"\n";
  const CutInterfaceConstructor* constructor = dynamic_cast<const CutInterfaceConstructor*>(geometryconstructor);
  assert(constructor);
  const alat::VectorOneVariable& phi = *constructor->getPhi();

  std::shared_ptr<const mesh::MeasureOfCell> measureofcell = std::dynamic_pointer_cast<const mesh::MeasureOfCell>(mesh->getGeometryObject(meshEnums::MeasureOfCell));
  assert(measureofcell);
  const arma::vec& moc = measureofcell->getMeasureOfCell();

  const EdgesAndCells& edgesandcells = mesh->getEdgesAndCells();
  const alat::armaimat& nodesofedges = edgesandcells.getNodesOfEdges();
  const alat::armaimat& edgesofcells = edgesandcells.getEdgesOfCells();
  const NodesAndNodesOfCells& nodesandnodesofcells = mesh->getNodesAndNodesOfCells();
  const arma::mat& nodes = nodesandnodesofcells.getNodes();
  const alat::armaimat& nodesofcells = nodesandnodesofcells.getNodesOfCells();
  assert(mesh->geometryObjectExists(meshEnums::Normals));
  std::shared_ptr<const mesh::Normals> normalsobj = std::dynamic_pointer_cast<const mesh::Normals>(mesh->getGeometryObject(meshEnums::Normals));
  assert(normalsobj);
  const arma::mat& normals = normalsobj->getNormals();
  const arma::fmat& sigmas = normalsobj->getSigma();
  const SidesAndCells& sidesandcells=mesh->getSidesAndCells();
  const alat::armaimat& nodesofsides = sidesandcells.getNodesOfSides();
  const alat::armaimat& sidesofcells = sidesandcells.getSidesOfCells();

  int nedgespercell = mesh->getNEdgesPerCell();
  int nnodespercell = mesh->getNNodesPerCell();
  int nedges = mesh->getNEdges();
  int nnodes = mesh->getNNodes();
  assert(nodesofedges.n_cols==nedges);
  int ncells = mesh->getNCells();
  _edgeiscut.set_size(nedges);
  _edgeiscut.fill(-1);
  _nodeiscut.set_size(nnodes);
  _nodeiscut.fill(-1);
  int ncutedges=0;
  for(int i=0;i<nedges;i++)
  {
    double phi0 = phi[nodesofedges(0,i)];
    double phi1 = phi[nodesofedges(1,i)];
    if(phi0*phi1 < 0)//-(fabs(phi0)+fabs(phi1)))
    {
      _edgeiscut[i] = ncutedges;
      ncutedges++;
    }
  }
  // std::cerr << "ncutedges="<<ncutedges<<"\n";
  _cutedges.set_size(ncutedges);
  _cutcoeff.set_size(ncutedges);
  for(int i=0;i<nedges;i++)
  {
    int iloc = _edgeiscut[i];
    if(iloc>=0)
    {
      _cutedges[iloc] = i;
      double phi0 = phi[nodesofedges(0,i)];
      double phi1 = phi[nodesofedges(1,i)];
      _cutcoeff[iloc] = phi0/(phi0-phi1);
    }
  }
  // std::cerr << "_edgeiscut="<<_edgeiscut.t();
  // std::cerr << "_cutedges="<<_cutedges.t();
  // std::cerr << "_cutcoeff="<<_cutcoeff.t();
  _celliscut.set_size(ncells);
  _celliscut.fill(-1);
  int ncutcells=0;
  for(int i=0;i<ncells;i++)
  {
    bool iscut=false;
    // std::cerr <<"iK="<<i << " edgesofcells=" << alat::armaivec(edgesofcells.col(i)).t();
    for(int ii=0;ii<nedgespercell;ii++)
    {
      if(_edgeiscut[edgesofcells(ii,i)]>=0)
      {
        iscut=true;
        break;
      }
    }
    if(iscut)
    {
      _celliscut[i] = ncutcells;
      ncutcells++;
    }
    else
    {
      if(phi[nodesofcells(0,i)]>0) {_celliscut[i] =-2;}
    }
  }
  // std::cerr << "ncutcells="<<ncutcells<<"\n";
  _cutcells.set_size(ncutcells);
  for(int i=0;i<ncells;i++)
  {
    int iloc = _celliscut[i];
    if(iloc>=0)
    {
      _cutcells[iloc] = i;
    }
  }
  _normalsofcutcells.set_size(3,ncutcells);
  _measuresofcutcells.set_size(2,ncutcells);
  _nodesofcutcellsisin.set_size(nnodespercell, ncutcells);
  _nodesofcutcellsisin.fill(0);
  // centers of cut cells
  _cinofcutcells.set_size(3, ncutcells);
  _cexofcutcells.set_size(3, ncutcells);
  _cofinofcutcells.set_size(nnodespercell, ncutcells);
  _cofexofcutcells.set_size(nnodespercell, ncutcells);

  for(int i=0;i<ncutcells;i++)
  {
    int iK = _cutcells[i];
    int iKcut = _celliscut[iK];
    assert(iKcut>=0);
    assert(i==iKcut);
    int countin=0,countex=0;
    for(int ii=0;ii<nnodespercell;ii++)
    {
      int iN = nodesofcells(ii,iK);
      _nodesofcutcellsisin(ii,i) = (phi[iN]<0.0);
      if(_nodesofcutcellsisin(ii,i)) {countin++;}
    }
    arma::mat innodes(3,countin), exnodes(3,nnodespercell-countin);
    arma::ivec innodesind(countin), exnodesind(nnodespercell-countin);
    countin=0;
    for(int ii=0;ii<nnodespercell;ii++)
    {
      int iN = nodesofcells(ii,iK);
      // _nodesofcutcellsisin(ii,i) = (phi[iN]<0.0);
      if(_nodesofcutcellsisin(ii,i))
      {
        innodes.col(countin) = nodes.col(iN);
        innodesind[countin] = iN;
        countin++;
      }
      else
      {
        exnodes.col(countex) = nodes.col(iN);
        exnodesind[countex] = iN;
        countex++;
      }
    }
    int cutedges=0;
    // std::cerr << "iK="<<iK<< " nodes=" << alat::armaivec(nodesofcells.col(iK)).t() << " edges=" << alat::armaivec(edgesofcells.col(iK)).t();
    for(int ii=0;ii<nedgespercell;ii++)
    {
      int iE = edgesofcells(ii,iK);
      int ilocedge = _edgeiscut[iE];
      if(ilocedge>=0)
      {
        cutedges++;
      }
    }
    arma::mat addnodes(3,cutedges);
    arma::ivec addnodesedges(cutedges);
    cutedges=0;
    for(int ii=0;ii<nedgespercell;ii++)
    {
      int iE = edgesofcells(ii,iK);
      // std::cerr << "iK="<<iK<<" iE="<<iE<< " is cut "<<  _edgeiscut[iE] <<" nodes-edge " << alat::armaivec(nodesofedges.col(iE)).t();
      int ilocedge = _edgeiscut[iE];
      if(ilocedge>=0)
      {
        double lambda = _cutcoeff[ilocedge];
        int iN0 = nodesofedges(0,iE);
        int iN1 = nodesofedges(1,iE);
        // std::cerr << "iN0="<<iN0<< " iN1="<<iN1<< " lambda="<<lambda<<"\n";
        addnodes.col(cutedges) = (1-lambda)*nodes.col(iN0)+lambda*nodes.col(iN1);
        addnodesedges[cutedges] = iE;
        // std::cerr << "iK="<<iK<<" iN0="<<iN0<< " iN1="<<iN1<< " addnodes.col(cutedges)="<<(addnodes.col(cutedges)).t();
        cutedges++;
      }
    }
    _computeMeasuresOfCutCell(iK, mesh, _measuresofcutcells.col(iKcut), _normalsofcutcells.col(iKcut), _cinofcutcells.col(iKcut), _cexofcutcells.col(iKcut), innodes, exnodes, addnodes, moc[iK], innodesind, exnodesind, addnodesedges);
    // std::cerr << "addnodes="<<addnodes;
    // std::cerr << " countin="<<countin<< " cutedges="<<cutedges<<"\n";

    double d = 1.0/moc[iK];
    d/=(double)mesh->getDimension();
    arma::vec xtest(3);
    for(int ii=0;ii<nnodespercell;ii++)
    {
      int iN = nodesofcells(ii,iK);
      int iS = sidesofcells(ii, iK);
      arma::subview_col<double> normal = normals.col(iS);
      double sigma = sigmas(ii,iK);
      xtest = _cinofcutcells.col(iKcut) - nodes.col(iN);
      _cofinofcutcells(ii,iKcut) =  1.0 - d*sigma*arma::dot(xtest,normal);
      xtest = _cexofcutcells.col(iKcut) - nodes.col(iN);
      _cofexofcutcells(ii,iKcut) =  1.0 - d*sigma*arma::dot(xtest,normal);
    }
    double testin=0.0, testex=0.0;
    for(int ii=0;ii<nnodespercell;ii++)
    {
      testin+=_cofinofcutcells(ii,iKcut);
      testex+=_cofexofcutcells(ii,iKcut);
    }
    if(fabs(testin-1.0)>1e-12 or fabs(testex-1.0)>1e-12)
    {
      std::cerr << "testin="<<testin<<"\n";
      std::cerr << "testex="<<testex<<"\n";
      assert(0);
    }
    // // test centers of cut cells
    // arma::vec coeff(nedgespercell, arma::fill::zeros), coeff2(nedgespercell, arma::fill::zeros);
    // alat::armaivec count(nnodespercell, arma::fill::zeros);
    // for(int ii=0;ii<nedgespercell;ii++)
    // {
    //   int ie = edgesandcells._edges_of_cells(ii,iK);
    //   int ile = _edgeiscut[ie];
    //   if(ile<0) continue;
    //   double cutcoeff = _cutcoeff[ile];
    //   int i0 = edgesandcells._localnodes_of_edges_in_cells(0, ii, iK);
    //   int i1 = edgesandcells._localnodes_of_edges_in_cells(1, ii, iK);
    //   int iN0 = nodesofcells(i0,iK);
    //   int iN1 = nodesofcells(i1,iK);
    //   arma::vec xc,x0, x1;
    //   x0 = nodes.col(iN0);
    //   x1 = nodes.col(iN1);
    //   xc = (1.0-cutcoeff) *x0  + cutcoeff*x1;
    //   // std::cerr << "iK="<<iK << " iN0,IN1="<<iN0 << ","<<iN1 << " cutcoeff="<<cutcoeff<<" x0="<<x0[0]<<" "<<x0[1] << " x1="<<x1[0]<<" "<<x1[1] << " xc="<<xc[0]<<" "<<xc[1] << "\n";
    //   coeff[i0] += (1.0-cutcoeff)/3.0;
    //   coeff[i1] += cutcoeff/3.0;
    //   count[i0]++;
    //   count[i1]++;
    // }
    // double K = moc[iK];
    // double Kin = _measuresofcutcells(0,iKcut);
    // double Kex = _measuresofcutcells(1,iKcut);
    // // std::cerr << "Kin="<<Kin << " Kex="<<Kex<<"\n";
    // if(mesh->getDimension()==2)
    // {
    //   for(int ii=0;ii<nnodespercell;ii++)
    //   {
    //     if(count[ii]!=2) continue;
    //     coeff[ii] += 1./3.;
    //     if(_nodesofcutcellsisin(ii,iKcut))
    //     {
    //       assert( arma::norm( _cofinofcutcells.col(iKcut) - coeff) < 1e-12);
    //       for(int jj=0;jj<3;jj++)
    //       {
    //         coeff2[jj] = (K/3.0 - Kin*_cofinofcutcells(jj,iKcut))/Kex;
    //       }
    //       assert( arma::norm( _cofexofcutcells.col(iKcut) - coeff2) < 1e-12);
    //     }
    //     else
    //     {
    //       assert( arma::norm( _cofexofcutcells.col(iKcut) - coeff) < 1e-12);
    //       for(int jj=0;jj<3;jj++)
    //       {
    //         coeff2[jj] = (K/3.0 - Kex*_cofexofcutcells(jj,iKcut))/Kin;
    //       }
    //       assert( arma::norm( _cofinofcutcells.col(iKcut) - coeff2) < 1e-12);
    //     }
    //   }
    // }
    // else
    // {
    //   assert(0);
    // }
    // std::cerr << "iK=" << iK<< " K=" << K << " Kin="<<Kin << " Kex="<<Kex<<"\n";
    // std::cerr << "coeffin="<<coeffin.t();
    // std::cerr << "coeffex="<<coeffex.t();
    // arma::vec xK(3, arma::fill::zeros), xin(3, arma::fill::zeros), xex(3, arma::fill::zeros);
    // xin.zeros(); xex.zeros();
    // for(int ii=0;ii<nnc;ii++)
    // {
    //   int iN = nodesofcells(ii,iK);
    //   xK += nodes.col(iN)/(double)nnc;
    //   xin += _cofinofcutcells(ii,iKcut) * nodes.col(iN);
    //   xex += _cofexofcutcells(ii,iKcut) * nodes.col(iN);
    //   // xin += coeffin[ii] * nodes.col(iN);
    //   // xex += coeffex[ii] * nodes.col(iN);
    // }
    // //test
    // // double pin = Kin/K;
    // // double pex = Kex/K;
    // // std::cerr << "iK="<<iK << " " << pin*coeffin+pex*coeffex;
    // arma::vec xtest = (Kin*xin+Kex*xex)/(Kin+Kex);
    // // std::cerr << "iK="<<iK << " xK="<<xK[0]<<" "<<xK[1] << " test="<<xtest[0]<<" "<<xtest[1] << "\n";
    // assert(arma::norm(xK-xtest)<1e-12);
    // _cinofcutcells.col(iKcut) = xin;
    // _cexofcutcells.col(iKcut) = xex;



    // // centers of cut cells
    // int nec = mesh->getNEdgesPerCell();
    // int nnc = mesh->getNNodesPerCell();
    // arma::vec coeff(nnc, arma::fill::zeros);
    // alat::armaivec count(nnc, arma::fill::zeros);
    // for(int ii=0;ii<nec;ii++)
    // {
    //   int ie = edgesandcells._edges_of_cells(ii,iK);
    //   int ile = _edgeiscut[ie];
    //   if(ile<0) continue;
    //   double cutcoeff = _cutcoeff[ile];
    //   int i0 = edgesandcells._localnodes_of_edges_in_cells(0, ii, iK);
    //   int i1 = edgesandcells._localnodes_of_edges_in_cells(1, ii, iK);
    //   int iN0 = nodesofcells(i0,iK);
    //   int iN1 = nodesofcells(i1,iK);
    //   arma::vec xc,x0, x1;
    //   x0 = nodes.col(iN0);
    //   x1 = nodes.col(iN1);
    //   xc = (1.0-cutcoeff) *x0  + cutcoeff*x1;
    //   // std::cerr << "iK="<<iK << " iN0,IN1="<<iN0 << ","<<iN1 << " cutcoeff="<<cutcoeff<<" x0="<<x0[0]<<" "<<x0[1] << " x1="<<x1[0]<<" "<<x1[1] << " xc="<<xc[0]<<" "<<xc[1] << "\n";
    //   coeff[i0] += (1.0-cutcoeff)/3.0;
    //   coeff[i1] += cutcoeff/3.0;
    //   count[i0]++;
    //   count[i1]++;
    // }
    // double K = moc[iK];
    // double Kin = _measuresofcutcells(0,iKcut);
    // double Kex = _measuresofcutcells(1,iKcut);
    // // std::cerr << "Kin="<<Kin << " Kex="<<Kex<<"\n";
    // if(mesh->getDimension()==2)
    // {
    //   for(int ii=0;ii<nnc;ii++)
    //   {
    //     if(count[ii]!=2) continue;
    //     coeff[ii] += 1./3.;
    //     if(_nodesofcutcellsisin(ii,iKcut))
    //     {
    //       _cofinofcutcells.col(iKcut) = coeff;
    //       for(int jj=0;jj<3;jj++)
    //       {
    //         _cofexofcutcells(jj,iKcut) = (K/3.0 - Kin*_cofinofcutcells(jj,iKcut))/Kex;
    //       }
    //     }
    //     else
    //     {
    //       _cofexofcutcells.col(iKcut) = coeff;
    //       for(int jj=0;jj<3;jj++)
    //       {
    //         _cofinofcutcells(jj,iKcut) = (K/3.0 - Kex*_cofexofcutcells(jj,iKcut))/Kin;
    //       }
    //     }
    //   }
    // }
    // else
    // {
    //   bool found=false;
    //   for(int ii=0;ii<nnc;ii++)
    //   {
    //     if(count[ii]==3)
    //     {
    //       found=true;
    //       coeff[ii] += 1./4.;
    //       if(_nodesofcutcellsisin(ii,iKcut))
    //       {
    //         _cofinofcutcells.col(iKcut) = coeff;
    //         for(int jj=0;jj<4;jj++)
    //         {
    //           _cofexofcutcells(jj,iKcut) = (K/4.0 - Kin*_cofinofcutcells(jj,iKcut))/Kex;
    //         }
    //       }
    //       else
    //       {
    //         _cofexofcutcells.col(iKcut) = coeff;
    //         for(int jj=0;jj<4;jj++)
    //         {
    //           _cofinofcutcells(jj,iKcut) = (K/4.0 - Kex*_cofexofcutcells(jj,iKcut))/Kin;
    //         }
    //       }
    //     }
    //   }
    //   assert(found);
    // }
    // // std::cerr << "iK=" << iK<< " K=" << K << " Kin="<<Kin << " Kex="<<Kex<<"\n";
    // // std::cerr << "coeffin="<<coeffin.t();
    // // std::cerr << "coeffex="<<coeffex.t();
    // arma::vec xK(3, arma::fill::zeros), xin(3, arma::fill::zeros), xex(3, arma::fill::zeros);
    // xin.zeros(); xex.zeros();
    // for(int ii=0;ii<nnc;ii++)
    // {
    //   int iN = nodesofcells(ii,iK);
    //   xK += nodes.col(iN)/(double)nnc;
    //   xin += _cofinofcutcells(ii,iKcut) * nodes.col(iN);
    //   xex += _cofexofcutcells(ii,iKcut) * nodes.col(iN);
    //   // xin += coeffin[ii] * nodes.col(iN);
    //   // xex += coeffex[ii] * nodes.col(iN);
    // }
    // //test
    // // double pin = Kin/K;
    // // double pex = Kex/K;
    // // std::cerr << "iK="<<iK << " " << pin*coeffin+pex*coeffex;
    // arma::vec xtest = (Kin*xin+Kex*xex)/(Kin+Kex);
    // // std::cerr << "iK="<<iK << " xK="<<xK[0]<<" "<<xK[1] << " test="<<xtest[0]<<" "<<xtest[1] << "\n";
    // assert(arma::norm(xK-xtest)<1e-12);
    // _cinofcutcells.col(iKcut) = xin;
    // _cexofcutcells.col(iKcut) = xex;

    // double d = 1.0/moc[iK];
    // if(mesh->getDimension()==2) d/=2.0;
    // else d/=6.0;
    // bool wrong=false;
    // for(int ii=0;ii<nnc;ii++)
    // {
    //   int iN = nodesofcells(ii,iK);
    //   int iS = sidesofcells(ii, iK);
    //   arma::subview_col<double> normal = normals.col(iS);
    //   double sigma = sigmas(ii,iK);
    //   xtest = xin - nodes.col(iN);
    //   double coefin =  1.0 - d*sigma*arma::dot(xtest,normal);
    //   if( fabs(_cofinofcutcells(ii,iKcut) - coefin) > 1e-12)
    //   {
    //     std::cerr << "_cofinofcutcells=" << _cofinofcutcells(ii,iKcut) << " coefin="<<coefin<<"\n";
    //     wrong=true;
    //   }
    // }
    // assert(not wrong);
  }
  // std::cerr << "_measuresofcutcells="<<_measuresofcutcells;
  // std::cerr << "_nodesofcutcellsisin="<<_nodesofcutcellsisin.t();
  // std::cerr << "_celliscut="<<_celliscut.t();
  // std::cerr << "_cutcells="<<_cutcells.t();






  alat::IntMap node2cutnodes;
  int count=0;
  for(int i=0;i<ncutcells;i++)
  {
    int iK = _cutcells[i];
    for(int ii=0;ii<nnodespercell;ii++)
    {
      int iN = nodesofcells(ii,iK);
      if(node2cutnodes.find(iN)==node2cutnodes.end())
      {
        node2cutnodes[iN] = count++;
      }
    }
  }
  _nodesofcutcell.set_size(nnodespercell, ncutcells);
  _cutnodes.set_size(count);
  _cutnodesisin.set_size(count);
  _cutnodesisin.zeros();
  for(alat::IntMap::const_iterator p=node2cutnodes.begin(); p!=node2cutnodes.end();p++)
  {
    _nodeiscut[p->first] = p->second;
    _cutnodes[p->second] = p->first;
  }
  for(int i=0;i<ncutcells;i++)
  {
    int iK = _cutcells[i];
    for(int ii=0;ii<nnodespercell;ii++)
    {
      int iN = nodesofcells(ii,iK);
      int iNcut = node2cutnodes[iN];
      assert(iNcut==_nodeiscut[iN]);
      if(_nodesofcutcellsisin(ii,i)) _cutnodesisin[iNcut] = 1;
      _nodesofcutcell(ii, i) = iNcut;
    }
  }

  // // centers of cut cells
  // _cinofcutcells.set_size(3, ncutcells);
  // _cexofcutcells.set_size(3, ncutcells);
  // _cofinofcutcells.set_size(nnodespercell, ncutcells);
  // _cofexofcutcells.set_size(nnodespercell, ncutcells);
  // for(int i=0;i<ncutcells;i++)
  // {
  //   int iK = _cutcells[i];
  //   int iKcut = _celliscut[iK];
  //   assert(iKcut>=0);
  //   assert(i==iKcut);
  //
  //   int nec = mesh->getNEdgesPerCell();
  //   int nnc = mesh->getNNodesPerCell();
  //   arma::vec coeff(nnc, arma::fill::zeros);
  //   alat::armaivec count(nnc, arma::fill::zeros);
  //   for(int ii=0;ii<nec;ii++)
  //   {
  //     int ie = edgesandcells._edges_of_cells(ii,iK);
  //     int ile = _edgeiscut[ie];
  //     if(ile<0) continue;
  //     double cutcoeff = _cutcoeff[ile];
  //     int i0 = edgesandcells._localnodes_of_edges_in_cells(0, ii, iK);
  //     int i1 = edgesandcells._localnodes_of_edges_in_cells(1, ii, iK);
  //     int iN0 = nodesofcells(i0,iK);
  //     int iN1 = nodesofcells(i1,iK);
  //     arma::vec xc,x0, x1;
  //     x0 = nodes.col(iN0);
  //     x1 = nodes.col(iN1);
  //     xc = (1.0-cutcoeff) *x0  + cutcoeff*x1;
  //     // std::cerr << "iK="<<iK << " iN0,IN1="<<iN0 << ","<<iN1 << " cutcoeff="<<cutcoeff<<" x0="<<x0[0]<<" "<<x0[1] << " x1="<<x1[0]<<" "<<x1[1] << " xc="<<xc[0]<<" "<<xc[1] << "\n";
  //     coeff[i0] += (1.0-cutcoeff)/3.0;
  //     coeff[i1] += cutcoeff/3.0;
  //     count[i0]++;
  //     count[i1]++;
  //   }
  //   double K = moc[iK];
  //   double Kin = _measuresofcutcells(0,iKcut);
  //   double Kex = _measuresofcutcells(1,iKcut);
  //   // std::cerr << "Kin="<<Kin << " Kex="<<Kex<<"\n";
  //   if(mesh->getDimension()==2)
  //   {
  //     for(int ii=0;ii<nnc;ii++)
  //     {
  //       if(count[ii]!=2) continue;
  //       coeff[ii] += 1./3.;
  //       if(_nodesofcutcellsisin(ii,iKcut))
  //       {
  //         _cofinofcutcells.col(iKcut) = coeff;
  //         for(int jj=0;jj<3;jj++)
  //         {
  //           _cofexofcutcells(jj,iKcut) = (K/3.0 - Kin*_cofinofcutcells(jj,iKcut))/Kex;
  //         }
  //       }
  //       else
  //       {
  //         _cofexofcutcells.col(iKcut) = coeff;
  //         for(int jj=0;jj<3;jj++)
  //         {
  //           _cofinofcutcells(jj,iKcut) = (K/3.0 - Kex*_cofexofcutcells(jj,iKcut))/Kin;
  //         }
  //       }
  //     }
  //   }
  //   else
  //   {
  //     bool found=false;
  //     for(int ii=0;ii<nnc;ii++)
  //     {
  //       if(count[ii]==3)
  //       {
  //         found=true;
  //         coeff[ii] += 1./4.;
  //         if(_nodesofcutcellsisin(ii,iKcut))
  //         {
  //           _cofinofcutcells.col(iKcut) = coeff;
  //           for(int jj=0;jj<4;jj++)
  //           {
  //             _cofexofcutcells(jj,iKcut) = (K/4.0 - Kin*_cofinofcutcells(jj,iKcut))/Kex;
  //           }
  //         }
  //         else
  //         {
  //           _cofexofcutcells.col(iKcut) = coeff;
  //           for(int jj=0;jj<4;jj++)
  //           {
  //             _cofinofcutcells(jj,iKcut) = (K/4.0 - Kex*_cofexofcutcells(jj,iKcut))/Kin;
  //           }
  //         }
  //       }
  //     }
  //     assert(found);
  //   }
  //
  //   // std::cerr << "iK=" << iK<< " K=" << K << " Kin="<<Kin << " Kex="<<Kex<<"\n";
  //   // std::cerr << "coeffin="<<coeffin.t();
  //   // std::cerr << "coeffex="<<coeffex.t();
  //   arma::vec xK(3, arma::fill::zeros), xin(3, arma::fill::zeros), xex(3, arma::fill::zeros);
  //   xin.zeros(); xex.zeros();
  //   for(int ii=0;ii<nnc;ii++)
  //   {
  //     int iN = nodesofcells(ii,iK);
  //     xK += nodes.col(iN)/(double)nnc;
  //     xin += _cofinofcutcells(ii,iKcut) * nodes.col(iN);
  //     xex += _cofexofcutcells(ii,iKcut) * nodes.col(iN);
  //     // xin += coeffin[ii] * nodes.col(iN);
  //     // xex += coeffex[ii] * nodes.col(iN);
  //   }
  //   //test
  //   // double pin = Kin/K;
  //   // double pex = Kex/K;
  //   // std::cerr << "iK="<<iK << " " << pin*coeffin+pex*coeffex;
  //   arma::vec xtest = (Kin*xin+Kex*xex)/(Kin+Kex);
  //   // std::cerr << "iK="<<iK << " xK="<<xK[0]<<" "<<xK[1] << " test="<<xtest[0]<<" "<<xtest[1] << "\n";
  //   assert(arma::norm(xK-xtest)<1e-12);
  //   _cinofcutcells.col(iKcut) = xin;
  //   _cexofcutcells.col(iKcut) = xex;
  //
  //   double d = 1.0/moc[iK];
  //   if(mesh->getDimension()==2) d/=2.0;
  //   else d/=6.0;
  //   bool wrong=false;
  //   for(int ii=0;ii<nnc;ii++)
  //   {
  //     int iN = nodesofcells(ii,iK);
  //     int iS = sidesofcells(ii, iK);
  //     arma::subview_col<double> normal = normals.col(iS);
  //     double sigma = sigmas(ii,iK);
  //     xtest = xin - nodes.col(iN);
  //     double coefin =  1.0 - d*sigma*arma::dot(xtest,normal);
  //     if( fabs(_cofinofcutcells(ii,iKcut) - coefin) > 1e-12)
  //     {
  //       std::cerr << "_cofinofcutcells=" << _cofinofcutcells(ii,iKcut) << " coefin="<<coefin<<"\n";
  //       wrong=true;
  //     }
  //   }
  //   assert(not wrong);
  // }

}

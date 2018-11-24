#include  "meshunitwithdata.hpp"
#include  <cassert>
#include  "p1cut.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
CutInterfaceConstructor::CutInterfaceConstructor(const alat::VectorOneVariable* vector) : mesh::CutInterfaceConstructor(), _vector(vector) {}
std::string CutInterfaceConstructor::getClassName() const {return "CutInterfaceConstructor";}
const alat::VectorOneVariable* CutInterfaceConstructor::getPhi() const {return _vector;}

/*--------------------------------------------------------------------------*/
MeshUnitWithData::~MeshUnitWithData() {}
MeshUnitWithData::MeshUnitWithData(std::string applicationname): solvers::MeshUnitWithData(),_applicationname(applicationname){}
MeshUnitWithData::MeshUnitWithData( const MeshUnitWithData& meshunitwithdata): solvers::MeshUnitWithData(meshunitwithdata)
{
  assert(0);
}
MeshUnitWithData& MeshUnitWithData::operator=( const MeshUnitWithData& meshunitwithdata)
{
  assert(0);
  solvers::MeshUnitWithData::operator=(meshunitwithdata);
  return *this;
}
std::string MeshUnitWithData::getClassName() const
{
  return "MeshUnitWithData";
}

/*--------------------------------------------------------------------------*/
void MeshUnitWithData::initDataVector(alat::VectorOneVariableInterface* u, int ivar, std::string varname) const
{
  if(varname=="Phi")
  {
    alat::VectorOneVariable* uv = dynamic_cast<alat::VectorOneVariable*>(u); assert(uv);
    uv->equal(&_phivector);
  }
  else
  {
    solvers::MeshUnitWithData::initDataVector(u, ivar, varname);
  }
}

/*--------------------------------------------------------------------------*/
void MeshUnitWithData::initMeshAndApplication(mesh::MeshUnitInterface* mesh, const solvers::SolverInterface* solver)
{
  solvers::MeshUnitWithData::initMeshAndApplication(mesh, solver);

  if(_mesh->geometryObjectExists(meshEnums::NodesCellsWeight))
  {
    assert(_mesh->geometryObjectExists(meshEnums::CutInterface));
    return;
  }

  // mesh->addGeometryObject(meshEnums::MeasureOfCell);
  // mesh->addGeometryObject(meshEnums::Normals);
  mesh->addGeometryObject(meshEnums::NodesCellsWeight);




  if(_applicationname=="linear_straight" or _applicationname=="quadratic_straight")
  {
    constructPhiByHand();
  }
  else
  {
    constructPhiBySplines(_applicationname);
  }
  // std::cerr << "_phivector = " << _phivector << "\n";



  CutInterfaceConstructor cutinterfaceconstructor(&_phivector);
  _mesh->addGeometryObject(meshEnums::CutInterface, &cutinterfaceconstructor);
  // std::cerr << "Mesh:\n" << _mesh->getInfo();

  assert(_mesh->geometryObjectExists(meshEnums::CutInterface));
  std::shared_ptr<const mesh::CutInterface> _cutinterface = std::dynamic_pointer_cast<const mesh::CutInterface>(_mesh->getGeometryObject(meshEnums::CutInterface));
  assert(_cutinterface);
  // std::cerr << "MeshUnitWithData::initMeshAndApplication() END\n";
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::constructPhiBySplines(std::string applicationname)
{
  int n = 8;
  alat::armamat x(3,n);
  if(applicationname=="cardioide")
  {
    double r = 0.5;
    for(int i=0;i<n;i++)
    {
      double phi = i*2*M_PI/n;
      double p = r * (1-cos(phi))*cos(phi)+0.5;
      double q = r * (1-cos(phi))*sin(phi);
      x(0, i) = -q;
      x(1, i) = p;
      x(2, i) = 0.0;
    }
  }
  else
  {
    double r = 0.5;
    for(int i=0;i<n;i++)
    {
      double phi = i*2*M_PI/n;
      x(0, i) = r * cos(phi);
      x(1, i) = r * sin(phi);
      x(2, i) = 0.0;
    }
  }
  _spline.generate(x);

  _phivector.set_size(_meshinfo->nnodes);
  for(int i=0;i<_meshinfo->nnodes;i++)
  {
    _phivector[i] = _spline.signeddistance(_meshinfo->nodes.col(i));
  }
  _spline.writeVtk("spline.vtk");
}



/*--------------------------------------------------------------------------*/
double MeshUnitWithData::_residualcorrectPhi(const solvers::FunctionInterface& phi, const alat::VectorOneVariable& phivector, int& ncutedges) const
{
  double res = 0.0;
  alat::armavec x0(3), x1(3), xbeta(3), uphi(1);
  ncutedges=0;
  for(int i=0;i<_meshinfo->nedges;i++)
  {
    int i0 = _meshinfo->nodes_of_edges(0,i);
    int i1 = _meshinfo->nodes_of_edges(1,i);
    double phi0 = phivector[i0];
    double phi1 = phivector[i1];
    if(phi0*phi1 >= 0) {continue;}
    ncutedges++;
    x0 = _meshinfo->nodes.col(i0);
    x1 = _meshinfo->nodes.col(i1);
    double beta = phi0/(phi0-phi1);
    xbeta = (1.0-beta)*x0 + beta*x1;
    phi(uphi, xbeta[0], xbeta[1], xbeta[2]);
    res += arma::dot(uphi,uphi);
  }
  return sqrt(res);
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::constructPhiByHand()
{
  // std::cerr << "MeshUnitWithData::constructPhi() BEGIN\n";
  const solvers::FunctionInterface& phi = _application->getDataFunction("Phi");
  _phivector.set_size(_meshinfo->nnodes);


  alat::armavec x0(3), x1(3), xbeta0(3), xbeta(3), xdiff(3), uphi(1), uphix(1), uphiy(1), uphiz(1);
  for(int i=0;i<_meshinfo->nnodes;i++)
  {
    phi(uphi, _meshinfo->nodes(0,i), _meshinfo->nodes(1,i), _meshinfo->nodes(2,i));
    _phivector[i] = uphi[0];
  }
  // std::cerr << "##### phi = " << _phivector << "\n";


  int ncutedges;
  _residualcorrectPhi(phi, _phivector, ncutedges);
  alat::armavec betagood(ncutedges);
  int count=0;
  for(int i=0;i<_meshinfo->nedges;i++)
  {
    int i0 = _meshinfo->nodes_of_edges(0,i);
    int i1 = _meshinfo->nodes_of_edges(1,i);
    double phi0 = _phivector[i0];
    double phi1 = _phivector[i1];
    if(phi0*phi1 >= 0) {continue;}
    ncutedges++;
    x0 = _meshinfo->nodes.col(i0);
    x1 = _meshinfo->nodes.col(i1);
    xdiff = x1-x0;
    double beta = phi0/(phi0-phi1);
    for(int k=0;k<10;k++)
    {
      xbeta = (1.0-beta)*x0 + beta*x1;
      phi(uphi, xbeta[0], xbeta[1], xbeta[2]);
      // std::cerr << " k="<<k << " beta="<<beta << " res=" << uphi[0] << "\n";
      if(fabs(uphi[0])<1e-12) {break;}
      phi.x(uphix, xbeta[0], xbeta[1], xbeta[2]);
      phi.y(uphiy, xbeta[0], xbeta[1], xbeta[2]);
      phi.z(uphiz, xbeta[0], xbeta[1], xbeta[2]);
      double grad = uphix[0]*xdiff[0]+uphiy[0]*xdiff[1];
      if(_mesh->getDimension()==3)
      {
        grad += uphiz[0]*xdiff[2];
      }
      double lambda = uphi[0]/grad;
      lambda = fmin(beta, fmax(beta-1.0, lambda) );
      beta -= lambda;
      // std::cerr << "beta="<<beta << " lambda="<<lambda<<"\n";
      assert(beta>=0.0);
      assert(beta<=1.0);
    }
    betagood[count] = beta;
    count++;
  }


  count=0;
  for(int i=0;i<_meshinfo->nedges;i++)
  {
    int i0 = _meshinfo->nodes_of_edges(0,i);
    int i1 = _meshinfo->nodes_of_edges(1,i);
    double phi0 = _phivector[i0];
    double phi1 = _phivector[i1];
    if(phi0*phi1 >= 0) {continue;}
    ncutedges++;
    x0 = _meshinfo->nodes.col(i0);
    x1 = _meshinfo->nodes.col(i1);
    xdiff = x1-x0;
    double beta0 = phi0/(phi0-phi1);
    double beta = betagood[count];
    xbeta0 = (1.0-beta0)*x0 + beta0*x1;
    xbeta = (1.0-beta)*x0 + beta*x1;
    // std::cerr << "[ "<<x0.t() << " , "<<x1.t()<<" ]" << " -> " << xbeta0.t() << " ? " << xbeta.t() << "\n";
    count++;
  }
  for(int k=0;k<1000;k++)
  {
    count=0;
    for(int i=0;i<_meshinfo->nedges;i++)
    {
      int i0 = _meshinfo->nodes_of_edges(0,i);
      int i1 = _meshinfo->nodes_of_edges(1,i);
      double phi0 = _phivector[i0];
      double phi1 = _phivector[i1];
      if(phi0*phi1 >= 0) {continue;}
      double beta = betagood[count];
      double lambda = (1.0-beta)*phi0 + beta*phi1;
      _phivector[i0] = phi0 - 0.9*lambda;
      _phivector[i1] = phi1 - 0.9*lambda;
      // std::cerr << "phi0="<<phi0<<" _phivector[i0]="<<_phivector[i0]<< "   phi1="<<phi1<<" _phivector[i1]="<<_phivector[i1]<<"\n";
      assert(_phivector[i0]*_phivector[i1]<0.0);
      assert(_phivector[i0]*phi0>0.0);
      assert(_phivector[i1]*phi1>0.0);
      count++;
    }
    double res = _residualcorrectPhi(phi, _phivector, ncutedges);
    // std::cerr << "k="<<k<<" res="<<res<<"\n";
    if(res<1e-10) {break;}
  }
  // std::cerr << "##### phi = " << _phivector << "\n";
  // assert(0);
  // std::cerr << "MeshUnitWithData::constructPhi() END\n";
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::_writeVtkCut(std::ofstream& file, const alat::armaivec& nodeintonode, const alat::IntMap& cellintocell, const alat::IntMap& nodetonodein, const alat::armavec& u)
{
  file<<"# vtk DataFile Version 4.0 "<<std::endl<<"output from SimFem"<<std::endl;
  file<<"ASCII"<<std::endl<<"DATASET UNSTRUCTURED_GRID"<<std::endl<<std::endl;

  int nn = nodeintonode.size();
  file<<"POINTS "<<nn<<" FLOAT"<<std::endl;
  for(int i = 0; i < nn; i++)
  {
    int iN = nodeintonode[i];
    file<<_meshinfo->nodes(0,iN)<<" "<<_meshinfo->nodes(1,iN)<<" "<<_meshinfo->nodes(2,iN)<<" "<<std::endl;
  }
  file<<std::endl;

  int ne = cellintocell.size();
  int nle = _meshinfo->nnodespercell;
  int lenght = ne*( nle+1 );

  file<<std::endl<<"CELLS "<<ne<<" "<<lenght<<std::endl;

  for(int ie = 0; ie < ne; ie++)
  {
    int iK = cellintocell[ie];
    file<<nle<<" ";
    for(int ii = 0; ii < nle; ii++)
    {
      int iN = _meshinfo->nodes_of_cells(ii, iK);
      file<< nodetonodein[iN] <<" ";
    }
    file<<std::endl;
  }
  file<<std::endl<<"CELL_TYPES "<<ne<<std::endl;
  for(int c = 0; c < ne; c++)
  {
    file<<5<<" ";
  }
  file<<std::endl << "POINT_DATA " << nn << "\n";
  file<<"SCALARS "<< "U"  << " double\n";
  file<<"LOOKUP_TABLE default\n";
  for(int i = 0; i < nn; i++)
  {
    file << u[i] << " ";
  }
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::writeVtkCut(std::string varname, std::string filename)
{
  alat::GhostVector gu("u");
  const alat::VectorAllVariables& u = *getVector(gu);
  assert(_mesh->geometryObjectExists(meshEnums::CutInterface));
  std::shared_ptr<const mesh::CutInterface> _cutinterface = std::dynamic_pointer_cast<const mesh::CutInterface>(_mesh->getGeometryObject(meshEnums::CutInterface));
  assert(_cutinterface);
  const alat::armaivec* cutcells = &_cutinterface->getCutCells();
  const alat::armaivec* celliscut = &_cutinterface->getCellIsCut();
  const alat::armaivec* cutnodes = &_cutinterface->getCutNodes();
  const alat::armaimat* nodesofcutcellsisin = &_cutinterface->getNodesOfCutCellsIsIn();
  const alat::armaimat* nodesofcutcells = &_cutinterface->getNodesOfCutCells();
  const alat::armaivec* nodeiscut = &_cutinterface->getNodeIsCut();
  const alat::armaivec* cutnodesisin = &_cutinterface->getCutNodeIsIn();

  int nnodes = _meshinfo->nnodes;
  alat::IntMap cellintocell, cellextocell;
  alat::IntMap nodetonodein, nodetonodeex;
  int countnodesin=0, countnodesex=0;
  int countcellsin=0, countcellsex=0;
  for(int iK=0;iK<_meshinfo->ncells;iK++)
  {
    int iKcut = (*celliscut)[iK];
    if(iKcut==-1)
    {
      cellintocell[countcellsin++]=iK;
      for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
      {
        int iN = _meshinfo->nodes_of_cells(ii,iK);
        if(nodetonodein.find(iN)==nodetonodein.end())
        {
          nodetonodein[iN] = countnodesin++;
        }
      }
    }
    else if(iKcut==-2)
    {
      cellextocell[countcellsex++]=iK;
      for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
      {
        int iN = _meshinfo->nodes_of_cells(ii,iK);
        if(nodetonodeex.find(iN)==nodetonodeex.end())
        {
          nodetonodeex[iN] = countnodesex++;
        }
      }
    }
    else if(iKcut>=0)
    {
      cellintocell[countcellsin++]=iK;
      cellextocell[countcellsex++]=iK;
      for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
      {
        int iN = _meshinfo->nodes_of_cells(ii,iK);
        if( (*nodesofcutcellsisin)(ii,iKcut) )
        {
          if(nodetonodein.find(iN)==nodetonodein.end())
          {
            nodetonodein[iN] = countnodesin++;
          }
          if(nodetonodeex.find(iN)==nodetonodeex.end())
          {
            nodetonodeex[iN] = countnodesex++;
          }
        }
        else
        {
          if(nodetonodein.find(iN)==nodetonodein.end())
          {
            nodetonodein[iN] = countnodesin++;
          }
          if(nodetonodeex.find(iN)==nodetonodeex.end())
          {
            nodetonodeex[iN] = countnodesex++;
          }
        }
      }
    }
  }
  alat::armaivec nodeintonode(countnodesin), nodeextonode(countnodesex);
  for(alat::IntMap::const_iterator p= nodetonodein.begin(); p!=nodetonodein.end();p++)
  {
    nodeintonode[p->second] = p->first;
  }
  for(alat::IntMap::const_iterator p= nodetonodeex.begin(); p!=nodetonodeex.end();p++)
  {
    nodeextonode[p->second] = p->first;
  }
  // std::cerr << "nodetonodein = " << nodetonodein << "\n";
  // std::cerr << "nodetonodeex = " << nodetonodeex << "\n";

  std::ofstream filein( (filename+"_in.vtk").c_str());
  std::ofstream fileex( (filename+"_ex.vtk").c_str());


  const alat::VectorOneVariable* u0 = dynamic_cast<const alat::VectorOneVariable*>(u.get(0));
  assert(u0);
  // std::cerr << "u0 = " << *u0 << "\n";
  // for(int iN=0;iN<nnodes;iN++)
  // {
  //   int iNcut = (*nodeiscut)[iN];
  //   if( iNcut>=0)
  //   {
  //     std::cerr << "iN iNcut = " << iN << " " << nnodes +iNcut << " : " << (*u0)[iN] << " " << (*u0)[nnodes + iNcut]<< "\n";
  //   }
  // }
  int ncomp = u0->ncomp();
  assert(ncomp==1);
  alat::armavec uin(nodeintonode.size()), uex(nodeextonode.size());

  for(int i=0; i< nodeintonode.size(); i++)
  {
    int iN = nodeintonode[i];
    int iNcut = (*nodeiscut)[iN];
    if( iNcut>=0  and (*cutnodesisin)[iNcut]==0 )
    {
      uin[i] = (*u0)[nnodes + iNcut];
      // std::cerr<<"x y uin " << iN << " " << nnodes +iNcut<< " " << i<< " : " << _meshinfo->nodes(0,iN)<<" "<<_meshinfo->nodes(1,iN)<<" -> "<<uin[i]<<"\n";
    }
    else
    {
      uin[i] = (*u0)[iN];
    }
  }
  for(int i=0; i< nodeextonode.size(); i++)
  {
    int iN = nodeextonode[i];
    int iNcut = (*nodeiscut)[iN];
    if( iNcut>=0 and (*cutnodesisin)[iNcut]==1)
    {
      uex[i] = (*u0)[nnodes + iNcut];
      // std::cerr<<"x y uex " << iN << " " << nnodes +iNcut<< " " << i<< " : " << _meshinfo->nodes(0,iN)<<" "<<_meshinfo->nodes(1,iN)<<" -> "<<uex[i]<<"\n";
    }
    else
    {
      uex[i] = (*u0)[iN];
    }
  }

  _writeVtkCut(filein, nodeintonode, cellintocell, nodetonodein, uin);
  _writeVtkCut(fileex, nodeextonode, cellextocell, nodetonodeex, uex);

  filein.close();
  fileex.close();
}
/*--------------------------------------------------------------------------*/
void MeshUnitWithData::writeVtkIIM(std::string varname, std::string filename)
{
  alat::GhostVector gu("u");
  const alat::VectorAllVariables& u = *getVector(gu);
  assert(_mesh->geometryObjectExists(meshEnums::CutInterface));
  std::shared_ptr<const mesh::CutInterface> _cutinterface = std::dynamic_pointer_cast<const mesh::CutInterface>(_mesh->getGeometryObject(meshEnums::CutInterface));
  assert(_cutinterface);
  const alat::armaivec* cutcells = &_cutinterface->getCutCells();
  const alat::armaivec* celliscut = &_cutinterface->getCellIsCut();
  const alat::armaivec* cutnodes = &_cutinterface->getCutNodes();
  const alat::armaimat* nodesofcutcellsisin = &_cutinterface->getNodesOfCutCellsIsIn();
  const alat::armaimat* nodesofcutcells = &_cutinterface->getNodesOfCutCells();
  const alat::armaivec* nodeiscut = &_cutinterface->getNodeIsCut();
  const alat::armaivec* cutnodesisin = &_cutinterface->getCutNodeIsIn();

  int nnodes = _meshinfo->nnodes;
  alat::IntMap cellintocell, cellextocell;
  alat::IntMap nodetonodein, nodetonodeex;
  int countnodesin=0, countnodesex=0;
  int countcellsin=0, countcellsex=0;
  for(int iK=0;iK<_meshinfo->ncells;iK++)
  {
    int iKcut = (*celliscut)[iK];
    if(iKcut==-1)
    {
      cellintocell[countcellsin++]=iK;
      for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
      {
        int iN = _meshinfo->nodes_of_cells(ii,iK);
        if(nodetonodein.find(iN)==nodetonodein.end())
        {
          nodetonodein[iN] = countnodesin++;
        }
      }
    }
    else if(iKcut==-2)
    {
      cellextocell[countcellsex++]=iK;
      for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
      {
        int iN = _meshinfo->nodes_of_cells(ii,iK);
        if(nodetonodeex.find(iN)==nodetonodeex.end())
        {
          nodetonodeex[iN] = countnodesex++;
        }
      }
    }
    // else if(iKcut>=0)
    // {
    //   cellintocell[countcellsin++]=iK;
    //   cellextocell[countcellsex++]=iK;
    //   for(int ii=0;ii<_meshinfo->nnodespercell;ii++)
    //   {
    //     int iN = _meshinfo->nodes_of_cells(ii,iK);
    //     if( (*nodesofcutcellsisin)(ii,iKcut) )
    //     {
    //       if(nodetonodein.find(iN)==nodetonodein.end())
    //       {
    //         nodetonodein[iN] = countnodesin++;
    //       }
    //       if(nodetonodeex.find(iN)==nodetonodeex.end())
    //       {
    //         nodetonodeex[iN] = countnodesex++;
    //       }
    //     }
    //     else
    //     {
    //       if(nodetonodein.find(iN)==nodetonodein.end())
    //       {
    //         nodetonodein[iN] = countnodesin++;
    //       }
    //       if(nodetonodeex.find(iN)==nodetonodeex.end())
    //       {
    //         nodetonodeex[iN] = countnodesex++;
    //       }
    //     }
    //   }
    // }
  }
  alat::armaivec nodeintonode(countnodesin), nodeextonode(countnodesex);
  for(alat::IntMap::const_iterator p= nodetonodein.begin(); p!=nodetonodein.end();p++)
  {
    nodeintonode[p->second] = p->first;
  }
  for(alat::IntMap::const_iterator p= nodetonodeex.begin(); p!=nodetonodeex.end();p++)
  {
    nodeextonode[p->second] = p->first;
  }
  // std::cerr << "nodetonodein = " << nodetonodein << "\n";
  // std::cerr << "nodetonodeex = " << nodetonodeex << "\n";

  std::ofstream filein( (filename+"_in.vtk").c_str());
  std::ofstream fileex( (filename+"_ex.vtk").c_str());


  const alat::VectorOneVariable* u0 = dynamic_cast<const alat::VectorOneVariable*>(u.get(0));
  assert(u0);
  // std::cerr << "u0 = " << *u0 << "\n";
  // for(int iN=0;iN<nnodes;iN++)
  // {
  //   int iNcut = (*nodeiscut)[iN];
  //   if( iNcut>=0)
  //   {
  //     std::cerr << "iN iNcut = " << iN << " " << nnodes +iNcut << " : " << (*u0)[iN] << " " << (*u0)[nnodes + iNcut]<< "\n";
  //   }
  // }
  int ncomp = u0->ncomp();
  assert(ncomp==1);
  alat::armavec uin(nodeintonode.size()), uex(nodeextonode.size());

  for(int i=0; i< nodeintonode.size(); i++)
  {
    int iN = nodeintonode[i];
    uin[i] = (*u0)[iN];
  }
  for(int i=0; i< nodeextonode.size(); i++)
  {
    int iN = nodeextonode[i];
    uex[i] = (*u0)[iN];
  }

  _writeVtkCut(filein, nodeintonode, cellintocell, nodetonodein, uin);
  _writeVtkCut(fileex, nodeextonode, cellextocell, nodetonodeex, uex);

  filein.close();
  fileex.close();
}

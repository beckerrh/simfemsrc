#include  "pdepartiim.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Alat/vectorallvariables.hpp"
#include  "Alat/matrixallvariables.hpp"

/*--------------------------------------------------------------------------*/
PdePartIIM::~PdePartIIM() {}
PdePartIIM::PdePartIIM(alat::StringList vars): PdePartCutInterface(vars){}
PdePartIIM::PdePartIIM( const PdePartIIM& pdepartwithfemtraditional): PdePartCutInterface(pdepartwithfemtraditional)
{
  assert(0);
}
PdePartIIM& PdePartIIM::operator=( const PdePartIIM& pdepartwithfemtraditional)
{
  assert(0);
  PdePartCutInterface::operator=(pdepartwithfemtraditional);
  return *this;
}
std::string PdePartIIM::getClassName() const
{
  return "PdePartIIM";
}
/*--------------------------------------------------------------------------*/
const solver_options::pdepart::opts PdePartIIM::setOptions()
{
  return solver_options::pdepart::cell+solver_options::pdepart::iside;
}

/*--------------------------------------------------------------------------*/
bool PdePartIIM::interiorsidecoupling(int iKin, int iKex) const
{
  if(p1iim or cr1iim)
  {
    if( (*_celliscut)(iKin)>=0 and (*_celliscut)(iKex)>=0) return true;
  }
  return false;
}

/*--------------------------------------------------------------------------*/
void PdePartIIM::setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters)
{
  PdePartCutInterface::setData(var2index, parameters);
  _gamma = parameters.doubles["gamma"];
  _gamma=0.0;
  p1iim = std::dynamic_pointer_cast<P1IIM>((*_fems)[_ivar]);
  p1iimex = std::dynamic_pointer_cast<P1IIM>((*_femsex)[_ivar]);
  cr1iim = std::dynamic_pointer_cast<CR1IIM>((*_fems)[_ivar]);
  cr1iimex = std::dynamic_pointer_cast<CR1IIM>((*_femsex)[_ivar]);
  if(p1iim)
  {
    p1iim->setCoefs(_kin, _kex);
    p1iimex->setCoefs(_kin, _kex);
  }
  if(cr1iim)
  {
    cr1iim->setCoefs(_kin, _kex);
    cr1iimex->setCoefs(_kin, _kex);
  }
  Linin.set_size(_nlocal,_nlocal);
  Linex.set_size(_nlocal,_nlocal);
  Lexin.set_size(_nlocal,_nlocal);
  Lexex.set_size(_nlocal,_nlocal);
  Last.set_size(_nlocal,_nlocal);
}
/*--------------------------------------------------------------------------*/
void PdePartIIM::computeRhsCellCut(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  arma::vec fin(_ncomp), fex(_ncomp);

  assert(_nlocal==3);
  floc[_ivar].zeros();
  int iKcut = (*_celliscut)[iK];
  double Kin = (*_measuresofcutcells)(0,iKcut);
  double Kex = (*_measuresofcutcells)(1,iKcut);
  const arma::subview_col<double> xin = (*_cinofcutcells).col(iKcut);
  const arma::subview_col<double> xex = (*_cexofcutcells).col(iKcut);
  const arma::subview_col<double> coeffin = (*_cofinofcutcells).col(iKcut);
  const arma::subview_col<double> coeffex = (*_cofexofcutcells).col(iKcut);
  _application->getRightHandSide(_ivar)(fin, xin[0], xin[1], xin[2]);
  _application->getRightHandSide(_ivar)(fex, xex[0], xex[1], xex[2]);
  if(p1iim)
  {
    int iKcut = (*_celliscut)[iK];
    double Kin = (*_measuresofcutcells)(0,iKcut);
    double Kex = (*_measuresofcutcells)(1,iKcut);
    const arma::subview_col<double> xin = (*_cinofcutcells).col(iKcut);
    const arma::subview_col<double> xex = (*_cexofcutcells).col(iKcut);
    const arma::subview_col<double> coeffin = (*_cofinofcutcells).col(iKcut);
    const arma::subview_col<double> coeffex = (*_cofexofcutcells).col(iKcut);

    _application->getRightHandSide(_ivar)(fin, xin[0], xin[1], xin[2]);
    _application->getRightHandSide(_ivar)(fex, xex[0], xex[1], xex[2]);
    for(int ii=0; ii<_nlocal;ii++)
    {
      int iN = _meshinfo->nodes_of_cells(ii,iK);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        for(int jj=0;jj<_nlocal;jj++)
        {
          floc[_ivar](icomp,ii) += p1iim->femcoefin(ii,jj)*coeffin[jj]*fin[icomp]*Kin;
          floc[_ivar](icomp,ii) += p1iim->femcoefex(ii,jj)*coeffex[jj]*fex[icomp]*Kex;
        }
      }
    }
  }
  else if(cr1iim)
  {
    int iKcut = (*_celliscut)[iK];
    double Kin = (*cr1iim->_measuresofcutcells)(0,iKcut);
    double Kex = (*cr1iim->_measuresofcutcells)(1,iKcut);
    const arma::subview_col<double> xin = (*_cinofcutcells).col(iKcut);
    const arma::subview_col<double> xex = (*_cexofcutcells).col(iKcut);
    const arma::subview_col<double> coeffin = (*_cofinofcutcells).col(iKcut);
    const arma::subview_col<double> coeffex = (*_cofexofcutcells).col(iKcut);

    // cr1iim->computeCutCentersOfGravity(iK);
    _application->getRightHandSide(_ivar)(fin, xin[0], xin[1], xin[2]);
    _application->getRightHandSide(_ivar)(fex, xex[0], xex[1], xex[2]);
    for(int ii=0; ii<_nlocal;ii++)
    {
      int iN = _meshinfo->nodes_of_cells(ii,iK);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        for(int jj=0;jj<_nlocal;jj++)
        {
          floc[_ivar](icomp,ii) += cr1iim->femcoefin(ii,jj)*coeffin[jj]*fin[icomp]*Kin;
          floc[_ivar](icomp,ii) += cr1iim->femcoefex(ii,jj)*coeffex[jj]*fex[icomp]*Kex;
        }
      }
    }
  }
  else
  {
    _error_string("computeRhsCellCut","unknown fem");
  }
}

/*--------------------------------------------------------------------------*/
void PdePartIIM::computeMatrixIIM(int iK)const
{
  if(p1iim) Last = p1iim->laplaceCut;
  else Last = cr1iim->laplaceCut;
}

/*--------------------------------------------------------------------------*/
void PdePartIIM::computeResidualCellCut(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  const arma::subview_row<double> uuloc = uloc[_ivar].row(0);
  arma::subview_row<double> ffloc = floc[_ivar].row(0);
  Last.zeros();
  computeMatrixIIM(iK);
  arma::vec fin = Last*uuloc.t();
  for(int ii=0; ii<_nlocal;ii++)
  {
    floc[_ivar](0,ii) += fin[ii];
  }
}
/*--------------------------------------------------------------------------*/
void PdePartIIM::computeMatrixCellCut(int iK, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const
{
  Last.zeros();
  computeMatrixIIM(iK);
  mat(_ivar,_ivar).set_size(_nlocal*_nlocal);
  mat_i(_ivar,_ivar).set_size(_nlocal*_nlocal);
  mat_j(_ivar,_ivar).set_size(_nlocal*_nlocal);
  mat(_ivar,_ivar).zeros();
  alat::armaivec indices;
  if(p1iim) p1iim->indicesOfCell(iK, indices);
  else cr1iim->indicesOfCell(iK, indices);
  assert(_ncomp==1);
  int count=0;
  for(int ii=0; ii<_nlocal;ii++)
  {
    for(int jj=0; jj<_nlocal;jj++)
    {
      mat(_ivar,_ivar)[count] += Last(ii,jj);
      mat_i(_ivar,_ivar)[count] = indices[ii];
      mat_j(_ivar,_ivar)[count] = indices[jj];
      count++;
    }
  }
}
/*--------------------------------------------------------------------------*/
void PdePartIIM::computeMatrixStabIIM(int iS, int iKin, int iKex)const
{
  const mesh::EdgesAndCells& edgesandcells = _mesh->getEdgesAndCells();
  int iN0 = _meshinfo->nodes_of_sides(0,iS);
  int iN1 = _meshinfo->nodes_of_sides(1,iS);
  int nec = _mesh->getNEdgesPerCell();
  int i0in, i1in, i0ex, i1ex;
  double cutcoefin, cutcoefex;
  for(int iie=0;iie<nec;iie++)
  {
    int ie = edgesandcells._edges_of_cells(iie,iKin);
    int iE0 = edgesandcells._nodes_of_edges(0,ie);
    int iE1 = edgesandcells._nodes_of_edges(1,ie);
    if( (iE0==iN0 and iE1==iN1) or (iE0==iN1 and iE1==iN0) )
    {
      i0in = edgesandcells._localnodes_of_edges_in_cells(0, iie, iKin);
      i1in = edgesandcells._localnodes_of_edges_in_cells(1, iie, iKin);
      int ile = (*_edgeiscut)[ie];
      if(ile>=0)
      {
        cutcoefin = (*_cutcoeff)[ile];
        break;
      }
    }
  }
  for(int iie=0;iie<nec;iie++)
  {
    int ie = edgesandcells._edges_of_cells(iie,iKex);
    int iE0 = edgesandcells._nodes_of_edges(0,ie);
    int iE1 = edgesandcells._nodes_of_edges(1,ie);
    if( (iE0==iN0 and iE1==iN1) or (iE0==iN1 and iE1==iN0) )
    {
      i0ex = edgesandcells._localnodes_of_edges_in_cells(0, iie, iKex);
      i1ex = edgesandcells._localnodes_of_edges_in_cells(1, iie, iKex);
      int ile = (*_edgeiscut)[ie];
      if(ile>=0)
      {
        cutcoefex = (*_cutcoeff)[ile];
        break;
      }
    }
  }
  assert(fabs(cutcoefin-cutcoefex)<1e-12);
  const arma::mat* femcoefin_in;
  const arma::mat* femcoefex_in;
  const arma::mat* femcoefin_ex;
  const arma::mat* femcoefex_ex;
  if(p1iim)
  {
    femcoefin_in = &p1iim->femcoefin;
    femcoefex_in = &p1iim->femcoefex;
    femcoefin_ex = &p1iimex->femcoefin;
    femcoefex_ex = &p1iimex->femcoefex;
  }
  else
  {
    femcoefin_in = &cr1iim->femcoefin;
    femcoefex_in = &cr1iim->femcoefex;
    femcoefin_ex = &cr1iimex->femcoefin;
    femcoefex_ex = &cr1iimex->femcoefex;
  }

  // for(int ii=0; ii<_nlocal;ii++)
  // {
  //   double din = (1.0-cutcoefin)*femcoefin_in(ii, i0in) + cutcoefin*femcoefin_in(ii, i1in);
  //   double dex = (1.0-cutcoefin)*femcoefex_in(ii, i0in) + cutcoefin*femcoefex_in(ii, i1in);
  //   // std::cerr << "IN din dex " << din << " " << dex << "\n";
  //   din = (1.0-cutcoefex)*femcoefin_ex(ii, i0ex) + cutcoefex*femcoefin_ex(ii, i1ex);
  //   dex = (1.0-cutcoefex)*femcoefex_ex(ii, i0ex) + cutcoefex*femcoefex_ex(ii, i1ex);
  //   // std::cerr << "EX din dex " << din << " " << dex << "\n";
  //   assert( fabs(din-dex)<1e-10);
  // }

  arma::vec coeffin(3), coeffex(3);
  coeffin.zeros(); coeffex.zeros();
  coeffin[i0in] = (1.0-cutcoefin);
  coeffin[i1in] = cutcoefin;
  coeffex[i0ex] = (1.0-cutcoefex);
  coeffex[i1ex] = cutcoefex;

  arma::vec lambdahatin_in = (*femcoefin_in)*coeffin;
  arma::vec lambdahatex_in = (*femcoefex_in)*coeffin;
  arma::vec lambdahatin_ex = (*femcoefin_ex)*coeffex;
  arma::vec lambdahatex_ex = (*femcoefex_ex)*coeffex;

  // double hinv = 1.0/sqrt(_meshinfo->measure_of_cells[iKin]);
  double gammain = _gamma*5.0*_localmodel->_kin;
  double gammaex = _gamma*5.0*_localmodel->_kex;
  Linin.zeros(); Linex.zeros(); Lexin.zeros(); Lexex.zeros();
  for(int ii=0; ii<_nlocal;ii++)
  {
    for(int jj=0; jj<_nlocal;jj++)
    {
      Linin(ii,jj) +=  gammain*lambdahatin_in[ii]*lambdahatin_in[jj];
      Linex(ii,jj) -=  gammain*lambdahatin_in[ii]*lambdahatin_ex[jj];
      Lexin(ii,jj) -=  gammain*lambdahatin_ex[ii]*lambdahatin_in[jj];
      Lexex(ii,jj) +=  gammain*lambdahatin_ex[ii]*lambdahatin_ex[jj];

      Linin(ii,jj) +=  gammaex*lambdahatex_in[ii]*lambdahatex_in[jj];
      Linex(ii,jj) -=  gammaex*lambdahatex_in[ii]*lambdahatex_ex[jj];
      Lexin(ii,jj) -=  gammaex*lambdahatex_ex[ii]*lambdahatex_in[jj];
      Lexex(ii,jj) +=  gammaex*lambdahatex_ex[ii]*lambdahatex_ex[jj];
    }
  }
}

/*--------------------------------------------------------------------------*/
void PdePartIIM::computeResidualInteriorSide(int iS, int iKin, int iKex, solvers::PdePartData::vec& flocin, solvers::PdePartData::vec& flocex, const solvers::PdePartData::vec& ulocin, const solvers::PdePartData::vec& ulocex)
{
  const arma::subview_row<double> uulocin = ulocin[_ivar].row(0);
  const arma::subview_row<double> uulocex = ulocex[_ivar].row(0);
  arma::subview_row<double> fflocin = flocin[_ivar].row(0);
  arma::subview_row<double> fflocex = flocex[_ivar].row(0);
  computeMatrixStabIIM(iS, iKin, iKex);
  arma::vec fin = Linin*uulocin.t()+Linex*uulocex.t();
  arma::vec fex = Lexin*uulocin.t()+Lexex*uulocex.t();
  // std::cerr << "Linin="<<Linin;
  // std::cerr << "uulocin="<<uulocin.t();
  // std::cerr << "uulocex="<<uulocex.t();
  // std::cerr << "fin="<<fin.t();
  // std::cerr << "fex="<<fex.t();
  for(int ii=0; ii<_nlocal;ii++)
  {
    flocin[_ivar](0,ii) += fin[ii];
    flocex[_ivar](0,ii) += fex[ii];
  }
}
/*--------------------------------------------------------------------------*/
void PdePartIIM::computeMatrixInteriorSide(int iS, int iKin, int iKex, solvers::PdePartData::mat& matinin, solvers::PdePartData::mat& matinex, solvers::PdePartData::mat& matexin, solvers::PdePartData::mat& matexex, solvers::PdePartData::imat& mat_i_in, solvers::PdePartData::imat& mat_j_in, solvers::PdePartData::imat& mat_i_ex, solvers::PdePartData::imat& mat_j_ex, const solvers::PdePartData::vec& ulocin, const solvers::PdePartData::vec& ulocex)const
{
  Linin.zeros(); Linex.zeros(); Lexin.zeros(); Lexex.zeros();
  computeMatrixStabIIM(iS, iKin, iKex);
  matinin(_ivar,_ivar).set_size(_nlocal*_nlocal);
  matinin(_ivar,_ivar).zeros();
  matinex(_ivar,_ivar).set_size(_nlocal*_nlocal);
  matinex(_ivar,_ivar).zeros();
  matexin(_ivar,_ivar).set_size(_nlocal*_nlocal);
  matexin(_ivar,_ivar).zeros();
  matexex(_ivar,_ivar).set_size(_nlocal*_nlocal);
  matexex(_ivar,_ivar).zeros();
  mat_i_in(_ivar,_ivar).set_size(_nlocal*_nlocal);
  mat_j_in(_ivar,_ivar).set_size(_nlocal*_nlocal);
  mat_i_ex(_ivar,_ivar).set_size(_nlocal*_nlocal);
  mat_j_ex(_ivar,_ivar).set_size(_nlocal*_nlocal);
  alat::armaivec indicesin, indicesex;
  (*_fems)[_ivar]->indicesOfCell(iKin, indicesin);
  (*_femsex)[_ivar]->indicesOfCell(iKex, indicesex);
  int count=0;
  for(int ii=0; ii<_nlocal;ii++)
  {
    for(int jj=0; jj<_nlocal;jj++)
    {
      matinin(_ivar,_ivar)[count] += Linin(ii,jj);
      matinex(_ivar,_ivar)[count] += Linex(ii,jj);
      matexin(_ivar,_ivar)[count] += Lexin(ii,jj);
      matexex(_ivar,_ivar)[count] += Lexex(ii,jj);
      mat_i_in(_ivar,_ivar)[count] = indicesin[ii];
      mat_j_in(_ivar,_ivar)[count] = indicesin[jj];
      mat_i_ex(_ivar,_ivar)[count] = indicesex[ii];
      mat_j_ex(_ivar,_ivar)[count] = indicesex[jj];
      count++;
    }
  }
}

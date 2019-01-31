#include  "pdepart.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Mesh/cutinterface.hpp"

/*--------------------------------------------------------------------------*/
PdePart::~PdePart() {}
PdePart::PdePart(alat::StringList vars): PdePartCutInterface(vars){}
PdePart::PdePart( const PdePart& pdepartwithfemtraditional): PdePartCutInterface(pdepartwithfemtraditional)
{
  assert(0);
}
PdePart& PdePart::operator=( const PdePart& pdepartwithfemtraditional)
{
  assert(0);
  PdePartCutInterface::operator=(pdepartwithfemtraditional);
  return *this;
}
std::string PdePart::getClassName() const
{
  return "PdePart";
}
/*--------------------------------------------------------------------------*/
const solver_options::pdepart::opts PdePart::setOptions()
{
  return solver_options::pdepart::cell;
}

/*--------------------------------------------------------------------------*/
void PdePart::setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters)
{
  PdePartCutInterface::setData(var2index, parameters);

  _gamma = parameters.doubles["gamma"];
  _method = parameters.strings["method"];
  p1cut = std::dynamic_pointer_cast<P1Cut>((*_fems)[_ivar]);
  p1hansbo = std::dynamic_pointer_cast<P1Hansbo>((*_fems)[_ivar]);
  if(p1cut)
  {
    p1cut->setCoefs(_kin, _kex);
  }
  if(p1hansbo)
  {
    p1hansbo->setCoefs(_kin, _kex);
  }
  if(_method=="nitsche" or _method=="strong")
  {
    assert(p1hansbo);
  }
  else
  {
    assert(p1cut);
  }
  Linin.set_size(_nlocal,_nlocal);
  Linex.set_size(_nlocal,_nlocal);
  Lexin.set_size(_nlocal,_nlocal);
  Lexex.set_size(_nlocal,_nlocal);
}
/*--------------------------------------------------------------------------*/
void PdePart::computeRhsCellCut(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const
{
  alat::armavec fin(_ncomp), fex(_ncomp);

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
  if(p1hansbo)
  {
    for(int ii=0; ii<_nlocal;ii++)
    {
      int iN = _meshinfo->nodes_of_cells(ii,iK);
      for(int icomp=0;icomp<_ncomp;icomp++)
      {
        floc[_ivar](icomp,p1hansbo->indin[ii]) += coeffin[ii]*fin[icomp]*Kin;
        floc[_ivar](icomp,p1hansbo->index[ii]) += coeffex[ii]*fex[icomp]*Kex;
      }
    }
  }
  else if(p1cut)
  {
    p1cut->computeBeta(iK);
    if(_method=="newnitsche")
    {
      for(int ii=0; ii<_nlocal;ii++)
      {
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          for(int jj=0;jj<_nlocal;jj++)
          {
            floc[_ivar](icomp,p1cut->indin[ii]) -= p1cut->P(ii,jj)*coeffin[jj]*fin[icomp]*Kin;
            floc[_ivar](icomp,p1cut->index[ii]) += p1cut->P(ii,jj)*coeffin[jj]*fin[icomp]*Kin;
            floc[_ivar](icomp,p1cut->index[ii]) -= p1cut->P(ii,jj)*coeffex[jj]*fex[icomp]*Kex;
            floc[_ivar](icomp,p1cut->indin[ii]) += p1cut->P(ii,jj)*coeffex[jj]*fex[icomp]*Kex;
            // floc[_ivar](icomp,p1cut->indin[ii]) -= Ain(ii,jj)*p1cut->coeffin[jj]*fin[icomp]*Kin;
            // floc[_ivar](icomp,p1cut->index[ii]) -= Aex(ii,jj)*p1cut->coeffex[jj]*fex[icomp]*Kex;
          }
        }
      }
    }
    else if(_method=="newnitsche2")
    {
      for(int ii=0; ii<_nlocal;ii++)
      {
        for(int icomp=0;icomp<_ncomp;icomp++)
        {
          for(int jj=0;jj<_nlocal;jj++)
          {
            floc[_ivar](icomp,p1cut->indin[ii]) -= p1cut->Pnc(ii,jj)*coeffin[jj]*fin[icomp]*Kin;
            floc[_ivar](icomp,p1cut->index[ii]) += p1cut->Pnc(ii,jj)*coeffin[jj]*fin[icomp]*Kin;
            floc[_ivar](icomp,p1cut->index[ii]) -= p1cut->Pnc(ii,jj)*coeffex[jj]*fex[icomp]*Kex;
            floc[_ivar](icomp,p1cut->indin[ii]) += p1cut->Pnc(ii,jj)*coeffex[jj]*fex[icomp]*Kex;
          }
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
void PdePart::computeResidualCellCut(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc) const
{
  const arma::subview_row<double> uuloc = uloc[_ivar].row(0);
  arma::subview_row<double> ffloc = floc[_ivar].row(0);

  Linin.zeros(); Linex.zeros(); Lexin.zeros(); Lexex.zeros();
  if(_method=="nitsche")
  {
    computeProjectionNitsche(iK);
  }
  else if(_method=="strong")
  {
    const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
    int iKcut = (*_celliscut)[iK];
    double moc=_meshinfo->measure_of_cells[iK];
    double Kin = (*_measuresofcutcells)(0,iKcut);
    double Kex = (*_measuresofcutcells)(1,iKcut);
    double diffin = _kin * Kin/moc;
    double diffex = _kex * Kex/moc;
    Linex.zeros(); Lexin.zeros();
    Linin = diffin*fem.laplace;
    Lexex = diffex*fem.laplace;
  }
  else if(_method=="nitschestabc")
  {
    computeProjectionNitscheStabC(iK);
  }
  else if(_method=="nitschestabnc")
  {
    computeProjectionNitscheStabNC(iK);
  }
  else if(_method=="newnitsche")
  {
    computeProjectionNewNitsche(iK);
  }
  else
  {
    computeProjectionNewNitsche2(iK);
  }
  const alat::armaivec* indin, *index;
  if(p1cut)
  {
    indin = &p1cut->indin;
    index = &p1cut->index;
  }
  else
  {
    indin = &p1hansbo->indin;
    index = &p1hansbo->index;
  }
  alat::armavec uin(_nlocal), uex(_nlocal);
  for(int ii=0; ii<_nlocal;ii++)
  {
    uin[ii] = uuloc[(*indin)[ii]];
    uex[ii] = uuloc[(*index)[ii]];
  }
  alat::armavec fin = Linin*uin + Linex*uex;
  alat::armavec fex = Lexin*uin + Lexex*uex;
  for(int ii=0; ii<_nlocal;ii++)
  {
    floc[_ivar](0,(*indin)[ii]) += fin[ii];
    floc[_ivar](0,(*index)[ii]) += fex[ii];
  }
}
/*--------------------------------------------------------------------------*/
void PdePart::computeMatrixCellCut(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const
{
  Linin.zeros(); Linex.zeros(); Lexin.zeros(); Lexex.zeros();
  if(_method=="nitsche")
  {
    computeProjectionNitsche(iK);
  }
  else if(_method=="strong")
  {
    const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
    int iKcut = (*_celliscut)[iK];
    double moc=_meshinfo->measure_of_cells[iK];
    double Kin = (*_measuresofcutcells)(0,iKcut);
    double Kex = (*_measuresofcutcells)(1,iKcut);
    double diffin = _kin * Kin/moc;
    double diffex = _kex * Kex/moc;
    Linex.zeros(); Lexin.zeros();
    Linin = diffin*fem.laplace;
    Lexex = diffex*fem.laplace;
  }
  else if(_method=="nitschestabc")
  {
    computeProjectionNitscheStabC(iK);
  }
  else if(_method=="nitschestabnc")
  {
    computeProjectionNitscheStabNC(iK);
  }
  else if(_method=="newnitsche")
  {
    computeProjectionNewNitsche(iK);
  }
  else
  {
    computeProjectionNewNitsche2(iK);
  }
  mat(_ivar,_ivar).set_size(4*_nlocal*_nlocal);
  mat_i(_ivar,_ivar).set_size(4*_nlocal*_nlocal);
  mat_j(_ivar,_ivar).set_size(4*_nlocal*_nlocal);
  mat(_ivar,_ivar).zeros();
  alat::armaivec indices;

  const alat::armaivec* indin, *index;
  if(p1cut)
  {
    indin = &p1cut->indin;
    index = &p1cut->index;
    p1cut->indicesOfCell(iK, indices);
  }
  else
  {
    indin = &p1hansbo->indin;
    index = &p1hansbo->index;
    p1hansbo->indicesOfCell(iK, indices);
  }

  assert(indices.size()==2*_nlocal);
  assert(_ncomp==1);
  int count=0;
  for(int ii=0; ii<_nlocal;ii++)
  {
    for(int jj=0; jj<_nlocal;jj++)
    {
      mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += Linin(ii,jj);
      mat_i(_ivar,_ivar)[count] = indices[(*indin)[ii]];
      mat_j(_ivar,_ivar)[count] = indices[(*indin)[jj]];
      count++;
      mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += Linex(ii,jj);
      mat_i(_ivar,_ivar)[count] = indices[(*indin)[ii]];
      mat_j(_ivar,_ivar)[count] = indices[(*index)[jj]];
      count++;
      mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += Lexin(ii,jj);
      mat_i(_ivar,_ivar)[count] = indices[(*index)[ii]];
      mat_j(_ivar,_ivar)[count] = indices[(*indin)[jj]];
      count++;
      mat(_ivar,_ivar)(icomp*_nlocal+ii, jcomp*_nlocal+jj) += Lexex(ii,jj);
      mat_i(_ivar,_ivar)[count] = indices[(*index)[ii]];
      mat_j(_ivar,_ivar)[count] = indices[(*index)[jj]];
      count++;
    }
  }
}

/*--------------------------------------------------------------------------*/
void PdePart::computeProjectionNitsche(int iK) const
{
  assert(_ncomp==1);
  int iKcut = (*_celliscut)[iK];
  double moc=_meshinfo->measure_of_cells[iK];
  double Kin = (*_measuresofcutcells)(0,iKcut);
  double Kex = (*_measuresofcutcells)(1,iKcut);
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  double diffin = _kin * Kin/moc;
  double diffex = _kex * Kex/moc;
  const arma::subview_col<double> normal = (*_normalsofcutcells).col(iKcut);
  double gammalength = arma::norm(normal);
  double gammagood = 2.0*_kin*_kex*gammalength/(_kin*Kex+_kex*Kin);
  double alphain = _kex*Kin/(_kin*Kex+_kex*Kin);
  gammagood *= 0.5*gammalength;
  gammagood*=_gamma;
  double alphaex = 1.0-alphain;
  int icomp=0;

  Linin.zeros(); Linex.zeros(); Lexin.zeros(); Lexex.zeros();
  Linin += diffin*fem.laplace;
  Lexex += diffex*fem.laplace;

  int nec = _mesh->getNEdgesPerCell();
  const mesh::EdgesAndCells& edgesandcells = _mesh->getEdgesAndCells();
  for(int ii=0;ii<nec;ii++)
  {
    int ie = edgesandcells._edges_of_cells(ii,iK);
    int ile = (*_edgeiscut)[ie];
    if(ile<0) continue;
    double cutcoeff = (*_cutcoeff)[ile];

    double d0 = (1.0-cutcoeff);
    double d1 = cutcoeff;

    double cici = d0*d0;
    double cice = d0*d1;
    double cece = d1*d1;

    int i0 = edgesandcells._localnodes_of_edges_in_cells(0, ii, iK);
    int i1 = edgesandcells._localnodes_of_edges_in_cells(1, ii, iK);

    Linin(i0, i0) += gammagood*cici;
    Linex(i0, i0) -= gammagood*cici;
    Lexin(i0, i0) -= gammagood*cici;
    Lexex(i0, i0) += gammagood*cici;
    Linin(i0, i1) += gammagood*cice;
    Linex(i0, i1) -= gammagood*cice;
    Lexin(i0, i1) -= gammagood*cice;
    Lexex(i0, i1) += gammagood*cice;
    Linin(i1, i0) += gammagood*cice;
    Linex(i1, i0) -= gammagood*cice;
    Lexin(i1, i0) -= gammagood*cice;
    Lexex(i1, i0) += gammagood*cice;
    Linin(i1, i1) += gammagood*cece;
    Linex(i1, i1) -= gammagood*cece;
    Lexin(i1, i1) -= gammagood*cece;
    Lexex(i1, i1) += gammagood*cece;


    for(int jj=0; jj<_nlocal;jj++)
    {
      double dphin = arma::dot(normal,fem.dphi.col(jj));
      Linin(i0, jj) -= 0.5*d0*_kin*alphain*dphin;
      Lexin(i0, jj) += 0.5*d0*_kin*alphain*dphin;
      Linin(i1, jj) -= 0.5*d1*_kin*alphain*dphin;
      Lexin(i1, jj) += 0.5*d1*_kin*alphain*dphin;

      Linex(i0, jj) -= 0.5*d0*_kex*alphaex*dphin;
      Lexex(i0, jj) += 0.5*d0*_kex*alphaex*dphin;
      Linex(i1, jj) -= 0.5*d1*_kex*alphaex*dphin;
      Lexex(i1, jj) += 0.5*d1*_kex*alphaex*dphin;

      Linin(jj, i0) -= 0.5*d0*_kin*alphain*dphin;
      Linex(jj, i0) += 0.5*d0*_kin*alphain*dphin;
      Linin(jj, i1) -= 0.5*d1*_kin*alphain*dphin;
      Linex(jj, i1) += 0.5*d1*_kin*alphain*dphin;

      Lexin(jj, i0) -= 0.5*d0*_kex*alphaex*dphin;
      Lexex(jj, i0) += 0.5*d0*_kex*alphaex*dphin;
      Lexin(jj, i1) -= 0.5*d1*_kex*alphaex*dphin;
      Lexex(jj, i1) += 0.5*d1*_kex*alphaex*dphin;
    }
  }
}

/*--------------------------------------------------------------------------*/
void PdePart::computeProjectionNitscheStabC(int iK) const
{
  assert(p1cut);
  p1cut->computeBeta(iK);
  assert(_ncomp==1);
  int iKcut = (*_celliscut)[iK];
  double moc=_meshinfo->measure_of_cells[iK];
  double Kin = (*_measuresofcutcells)(0,iKcut);
  double Kex = (*_measuresofcutcells)(1,iKcut);
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  double diffin = _kin * Kin/moc;
  double diffex = _kex * Kex/moc;
  alat::armaivec indices;
  p1cut->indicesOfCell(iK, indices);
  assert(indices.size()==2*_nlocal);
  const arma::subview_col<double> normal = (*_normalsofcutcells).col(iKcut);
  double gammalength = arma::norm(normal);
  double gammagood = 2.0*_kin*_kex*gammalength/(_kin*Kex+_kex*Kin);
  double alphain = _kex*Kin/(_kin*Kex+_kex*Kin);
  gammagood *= 0.5*gammalength;
  gammagood*=_gamma;
  double alphaex = 1.0-alphain;
  int nec = _mesh->getNEdgesPerCell();
  const mesh::EdgesAndCells& edgesandcells = _mesh->getEdgesAndCells();
  int icomp=0;

  Linin.zeros(); Linex.zeros(); Lexin.zeros(); Lexex.zeros();
  Linin += diffin*fem.laplace;
  Lexex += diffex*fem.laplace;

  double d = _gamma*(diffin+diffex);
  Linin += d*p1cut->P*fem.laplace*p1cut->P;
  Linex -= d*p1cut->P*fem.laplace*p1cut->P;
  Lexin -= d*p1cut->P*fem.laplace*p1cut->P;
  Lexex += d*p1cut->P*fem.laplace*p1cut->P;

  for(int ii=0;ii<nec;ii++)
  {
    int ie = edgesandcells._edges_of_cells(ii,iK);
    int ile = (*_edgeiscut)[ie];
    if(ile<0) continue;
    double cutcoeff = (*_cutcoeff)[ile];

    double d0 = (1.0-cutcoeff);
    double d1 = cutcoeff;

    int i0 = edgesandcells._localnodes_of_edges_in_cells(0, ii, iK);
    int i1 = edgesandcells._localnodes_of_edges_in_cells(1, ii, iK);
    for(int jj=0; jj<_nlocal;jj++)
    {
      double dphin = arma::dot(normal,fem.dphi.col(jj));
      Linin(i0, jj) -= 0.5*d0*_kin*alphain*dphin;
      Lexin(i0, jj) += 0.5*d0*_kin*alphain*dphin;
      Linin(i1, jj) -= 0.5*d1*_kin*alphain*dphin;
      Lexin(i1, jj) += 0.5*d1*_kin*alphain*dphin;

      Linin(jj, i0) -= 0.5*d0*_kin*alphain*dphin;
      Linex(jj, i0) += 0.5*d0*_kin*alphain*dphin;
      Linin(jj, i1) -= 0.5*d1*_kin*alphain*dphin;
      Linex(jj, i1) += 0.5*d1*_kin*alphain*dphin;

      Linex(i0, jj) -= 0.5*d0*_kex*alphaex*dphin;
      Lexex(i0, jj) += 0.5*d0*_kex*alphaex*dphin;
      Linex(i1, jj) -= 0.5*d1*_kex*alphaex*dphin;
      Lexex(i1, jj) += 0.5*d1*_kex*alphaex*dphin;

      Lexin(jj, i0) -= 0.5*d0*_kex*alphaex*dphin;
      Lexex(jj, i0) += 0.5*d0*_kex*alphaex*dphin;
      Lexin(jj, i1) -= 0.5*d1*_kex*alphaex*dphin;
      Lexex(jj, i1) += 0.5*d1*_kex*alphaex*dphin;
    }
  }
}

/*--------------------------------------------------------------------------*/
void PdePart::computeProjectionNitscheStabNC(int iK) const
{
  assert(p1cut);
  p1cut->computeBeta(iK);
  assert(_ncomp==1);
  int iKcut = (*_celliscut)[iK];
  double moc=_meshinfo->measure_of_cells[iK];
  double Kin = (*_measuresofcutcells)(0,iKcut);
  double Kex = (*_measuresofcutcells)(1,iKcut);
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  double diffin = _kin * Kin/moc;
  double diffex = _kex * Kex/moc;
  alat::armaivec indices;
  p1cut->indicesOfCell(iK, indices);
  assert(indices.size()==2*_nlocal);
  const arma::subview_col<double> normal = (*_normalsofcutcells).col(iKcut);
  double gammalength = arma::norm(normal);
  double gammagood = 2.0*_kin*_kex*gammalength/(_kin*Kex+_kex*Kin);
  double alphain = _kex*Kin/(_kin*Kex+_kex*Kin);
  gammagood *= 0.5*gammalength;
  gammagood*=_gamma;
  // alphain=0.0;
  double alphaex = 1.0-alphain;
  const mesh::EdgesAndCells& edgesandcells = _mesh->getEdgesAndCells();
  int icomp=0;

  Linin.zeros(); Linex.zeros(); Lexin.zeros(); Lexex.zeros();
  Linin += diffin*fem.laplace;
  Lexex += diffex*fem.laplace;

  double d = _gamma*(diffin+diffex);
  Linin += d*p1cut->Pnc*fem.laplace*p1cut->Pnc;
  Linex -= d*p1cut->Pnc*fem.laplace*p1cut->Pnc;
  Lexin -= d*p1cut->Pnc*fem.laplace*p1cut->Pnc;
  Lexex += d*p1cut->Pnc*fem.laplace*p1cut->Pnc;

  int nec = _mesh->getNEdgesPerCell();
  for(int ii=0;ii<nec;ii++)
  {
    int ie = edgesandcells._edges_of_cells(ii,iK);
    int ile = (*_edgeiscut)[ie];
    if(ile<0) continue;
    double cutcoeff = (*_cutcoeff)[ile];

    double d0 = (1.0-cutcoeff);
    double d1 = cutcoeff;

    int i0 = edgesandcells._localnodes_of_edges_in_cells(0, ii, iK);
    int i1 = edgesandcells._localnodes_of_edges_in_cells(1, ii, iK);

    for(int jj=0; jj<_nlocal;jj++)
    {
      double dphin = arma::dot(normal,fem.dphi.col(jj));
      Linin(i0, jj) -= 0.5*d0*_kin*alphain*dphin;
      Lexin(i0, jj) += 0.5*d0*_kin*alphain*dphin;
      Linin(i1, jj) -= 0.5*d1*_kin*alphain*dphin;
      Lexin(i1, jj) += 0.5*d1*_kin*alphain*dphin;

      Linin(jj, i0) -= 0.5*d0*_kin*alphain*dphin;
      Linex(jj, i0) += 0.5*d0*_kin*alphain*dphin;
      Linin(jj, i1) -= 0.5*d1*_kin*alphain*dphin;
      Linex(jj, i1) += 0.5*d1*_kin*alphain*dphin;

      Linex(i0, jj) -= 0.5*d0*_kex*alphaex*dphin;
      Lexex(i0, jj) += 0.5*d0*_kex*alphaex*dphin;
      Linex(i1, jj) -= 0.5*d1*_kex*alphaex*dphin;
      Lexex(i1, jj) += 0.5*d1*_kex*alphaex*dphin;

      Lexin(jj, i0) -= 0.5*d0*_kex*alphaex*dphin;
      Lexex(jj, i0) += 0.5*d0*_kex*alphaex*dphin;
      Lexin(jj, i1) -= 0.5*d1*_kex*alphaex*dphin;
      Lexex(jj, i1) += 0.5*d1*_kex*alphaex*dphin;
    }
  }
}


/*--------------------------------------------------------------------------*/
void PdePart::computeProjectionNewNitsche(int iK) const
{
  assert(p1cut);
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  int iKcut = (*_celliscut)[iK];
  double moc=_meshinfo->measure_of_cells[iK];
  double Kin = (*_measuresofcutcells)(0,iKcut);
  double Kex = (*_measuresofcutcells)(1,iKcut);
  double diffin = _kin * Kin/moc;
  double diffex = _kex * Kex/moc;

  p1cut->computeBeta(iK);

  // Linin = diffin*(p1cut->Iin + p1cut->Q*p1cut->Iex)*fem.laplace;
  // Lexex = diffex*(p1cut->Iex + p1cut->Q*p1cut->Iin)*fem.laplace;
  // Lexin = diffin*p1cut->P*p1cut->Iex*fem.laplace;
  // Linex = diffex*p1cut->P*p1cut->Iin*fem.laplace;

  Linin.zeros();  Linex.zeros(); Lexin.zeros();  Lexex.zeros();

    // Linin = diffin*fem.laplace;
    // Lexex = diffex*fem.laplace;

    Linin += diffin*( p1cut->Q)*fem.laplace*(p1cut->Q);
    Linex += diffin*( p1cut->Q)*fem.laplace*(p1cut->P);
    Lexin += diffin*( p1cut->P)*fem.laplace*(p1cut->Q);
    Lexex += diffin*( p1cut->P)*fem.laplace*(p1cut->P);

    Lexex += diffex*( p1cut->Q)*fem.laplace*(p1cut->Q);
    Lexin += diffex*( p1cut->Q)*fem.laplace*(p1cut->P);
    Linex += diffex*( p1cut->P)*fem.laplace*(p1cut->Q);
    Linin += diffex*( p1cut->P)*fem.laplace*(p1cut->P);


    double d = 8.0*_gamma*(diffin+diffex);
    // double d = (diffin+diffex);
    Linin += d*(p1cut->P*fem.laplace*p1cut->P);
    Linex -= d*(p1cut->P*fem.laplace*p1cut->P);
    Lexin -= d*(p1cut->P*fem.laplace*p1cut->P);
    Lexex += d*(p1cut->P*fem.laplace*p1cut->P);





  // Linin -= diffin*p1cut->Pnc*p1cut->Iex*fem.laplace;
  // Lexin += diffin*p1cut->Pnc*p1cut->Iex*fem.laplace;
  // Linex += diffex*p1cut->Pnc*p1cut->Iin*fem.laplace;
  // Lexex -= diffex*p1cut->Pnc*p1cut->Iin*fem.laplace;

  // Linin -= diffin*p1cut->P*p1cut->Iex*fem.laplace;
  // Lexin += diffin*p1cut->P*p1cut->Iex*fem.laplace;
  // Linex += diffex*p1cut->P*p1cut->Iin*fem.laplace;
  // Lexex -= diffex*p1cut->P*p1cut->Iin*fem.laplace;
  //
  // double d = 2.0*_gamma*(diffin+diffex);
  // Linin += d*(p1cut->P*fem.laplace*p1cut->P);
  // Linex -= d*(p1cut->P*fem.laplace*p1cut->P);
  // Lexin -= d*(p1cut->P*fem.laplace*p1cut->P);
  // Lexex += d*(p1cut->P*fem.laplace*p1cut->P);

  // double din = 2.0*_gamma*diffin;
  // Linin += din*(p1cut->P*p1cut->Iex*fem.laplace*p1cut->Iex*p1cut->P);
  // Linex -= din*(p1cut->P*p1cut->Iex*fem.laplace*p1cut->Iex*p1cut->P);
  // Lexin -= din*(p1cut->P*p1cut->Iex*fem.laplace*p1cut->Iex*p1cut->P);
  // Lexex += din*(p1cut->P*p1cut->Iex*fem.laplace*p1cut->Iex*p1cut->P);
  // double dex = 2.0*_gamma*diffex;
  // Linin += dex*(p1cut->P*p1cut->Iin*fem.laplace*p1cut->Iin*p1cut->P);
  // Linex -= dex*(p1cut->P*p1cut->Iin*fem.laplace*p1cut->Iin*p1cut->P);
  // Lexin -= dex*(p1cut->P*p1cut->Iin*fem.laplace*p1cut->Iin*p1cut->P);
  // Lexex += dex*(p1cut->P*p1cut->Iin*fem.laplace*p1cut->Iin*p1cut->P);

  // Linin += diffin*(p1cut->Iin*fem.laplace);
  // Lexin += diffin*(p1cut->Iex*fem.laplace);
  // Linex += diffex*(p1cut->Iin*fem.laplace);
  // Lexex += diffex*(p1cut->Iex*fem.laplace);



  // // enforce la continuite
  // Linin += diffin*p1cut->Iex*fem.laplace*p1cut->Iex;
  // Linex -= diffin*p1cut->Iex*fem.laplace*p1cut->Iex;
  // // Lexin -= diffin*p1cut->Iex*fem.laplace*p1cut->Iex;
  // // Lexex += diffin*p1cut->Iex*fem.laplace*p1cut->Iex;
  // // Linin += diffex*p1cut->Iin*fem.laplace*p1cut->Iin;
  // // Linex -= diffex*p1cut->Iin*fem.laplace*p1cut->Iin;
  // Lexin -= diffex*p1cut->Iin*fem.laplace*p1cut->Iin;
  // Lexex += diffex*p1cut->Iin*fem.laplace*p1cut->Iin;


  // const arma::subview_col<double> normal = (*_normalsofcutcells).col(iKcut);
  // double gammalength = arma::norm(normal);
  // double gammagood = 2.0*_kin*_kex*gammalength/(_kin*Kex+_kex*Kin);
  // double alphain = _kex*Kin/(_kin*Kex+_kex*Kin);
  // gammagood *= 0.5*gammalength;
  // gammagood*=_gamma;
  // int nec = _mesh->getNEdgesPerCell();
  // const mesh::EdgesAndCells& edgesandcells = _mesh->getEdgesAndCells();
  // for(int ii=0;ii<nec;ii++)
  // {
  //   int ie = edgesandcells._edges_of_cells(ii,iK);
  //   int ile = (*p1cut->_edgeiscut)[ie];
  //   if(ile<0) continue;
  //   double cutcoeff = (*_cutcoeff)[ile];
  //
  //   double d0 = (1.0-cutcoeff);
  //   double d1 = cutcoeff;
  //
  //   double cici = d0*d0;
  //   double cice = d0*d1;
  //   double cece = d1*d1;
  //
  //   int i0 = edgesandcells._localnodes_of_edges_in_cells(0, ii, iK);
  //   int i1 = edgesandcells._localnodes_of_edges_in_cells(1, ii, iK);
  //
  //   Linin(i0, i0) += gammagood*cici;
  //   Linex(i0, i0) -= gammagood*cici;
  //   Lexin(i0, i0) -= gammagood*cici;
  //   Lexex(i0, i0) += gammagood*cici;
  //   Linin(i0, i1) += gammagood*cice;
  //   Linex(i0, i1) -= gammagood*cice;
  //   Lexin(i0, i1) -= gammagood*cice;
  //   Lexex(i0, i1) += gammagood*cice;
  //   Linin(i1, i0) += gammagood*cice;
  //   Linex(i1, i0) -= gammagood*cice;
  //   Lexin(i1, i0) -= gammagood*cice;
  //   Lexex(i1, i0) += gammagood*cice;
  //   Linin(i1, i1) += gammagood*cece;
  //   Linex(i1, i1) -= gammagood*cece;
  //   Lexin(i1, i1) -= gammagood*cece;
  //   Lexex(i1, i1) += gammagood*cece;
  // }


  // Linin -= 0.5*diffin*p1cut->Pnc*p1cut->Iex*fem.laplace;
  // Lexin += 0.5*diffin*p1cut->Pnc*p1cut->Iex*fem.laplace;
  //
  // Linex += 0.5*diffex*p1cut->Pnc*p1cut->Iin*fem.laplace;
  // Lexex -= 0.5*diffex*p1cut->Pnc*p1cut->Iin*fem.laplace;



  // Linin -= diffin*p1cut->Pnc*p1cut->Iex*fem.laplace;
  // Lexin += diffin*p1cut->Pnc*p1cut->Iex*fem.laplace;
  // Linex += diffex*p1cut->Pnc*p1cut->Iin*fem.laplace;
  // Lexex -= diffex*p1cut->Pnc*p1cut->Iin*fem.laplace;
  //
  // Linin -= diffin*fem.laplace*p1cut->Iex*p1cut->Pnc;
  // Linex += diffin*fem.laplace*p1cut->Iex*p1cut->Pnc;
  // Lexin += diffex*fem.laplace*p1cut->Iin*p1cut->Pnc;
  // Lexex -= diffex*fem.laplace*p1cut->Iin*p1cut->Pnc;

  // double d = 2.0*_gamma*(diffin+diffex);
  // Linin += d*(p1cut->P*fem.laplace*p1cut->P);
  // Linex -= d*(p1cut->P*fem.laplace*p1cut->P);
  // Lexin -= d*(p1cut->P*fem.laplace*p1cut->P);
  // Lexex += d*(p1cut->P*fem.laplace*p1cut->P);
}

/*--------------------------------------------------------------------------*/
void PdePart::computeProjectionNewNitsche2(int iK) const
{
  assert(p1cut);
  const solvers::FemData& fem = (*_fems)[_ivar]->getFemdata();
  int iKcut = (*_celliscut)[iK];
  double moc=_meshinfo->measure_of_cells[iK];
  double Kin = (*_measuresofcutcells)(0,iKcut);
  double Kex = (*_measuresofcutcells)(1,iKcut);
  double diffin = _kin * Kin/moc;
  double diffex = _kex * Kex/moc;

  p1cut->computeBeta(iK);
  Linin.zeros();  Linex.zeros(); Lexin.zeros();  Lexex.zeros();

  // Linin = diffin*fem.laplace;
  // Lexex = diffex*fem.laplace;

  Linin += diffin*( p1cut->Qnc)*fem.laplace*(p1cut->Qnc);
  Linex += diffin*( p1cut->Qnc)*fem.laplace*(p1cut->Pnc);
  Lexin += diffin*( p1cut->Pnc)*fem.laplace*(p1cut->Qnc);
  Lexex += diffin*( p1cut->Pnc)*fem.laplace*(p1cut->Pnc);

  Lexex += diffex*( p1cut->Qnc)*fem.laplace*(p1cut->Qnc);
  Lexin += diffex*( p1cut->Qnc)*fem.laplace*(p1cut->Pnc);
  Linex += diffex*( p1cut->Pnc)*fem.laplace*(p1cut->Qnc);
  Linin += diffex*( p1cut->Pnc)*fem.laplace*(p1cut->Pnc);


  double d = 8.0*_gamma*(diffin+diffex);
  // double d = (diffin+diffex);
  Linin += d*(p1cut->P*fem.laplace*p1cut->P);
  Linex -= d*(p1cut->P*fem.laplace*p1cut->P);
  Lexin -= d*(p1cut->P*fem.laplace*p1cut->P);
  Lexex += d*(p1cut->P*fem.laplace*p1cut->P);
}

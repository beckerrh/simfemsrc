#include  "pdepartmultiplier.hpp"
#include  "Mesh/edgesandcells.hpp"
#include  "Mesh/cutinterface.hpp"
#include  "Alat/sparsitypatternsoft.hpp"
#include  "Alat/vectorallvariables.hpp"
#include  "Alat/matrixallvariables.hpp"

/*--------------------------------------------------------------------------*/
PdePartMultiplier::~PdePartMultiplier() {}
PdePartMultiplier::PdePartMultiplier(alat::StringList vars): PdePartCutInterface(vars){}
PdePartMultiplier::PdePartMultiplier( const PdePartMultiplier& pdepartwithfemtraditional): PdePartCutInterface(pdepartwithfemtraditional)
{
  assert(0);
}
PdePartMultiplier& PdePartMultiplier::operator=( const PdePartMultiplier& pdepartwithfemtraditional)
{
  assert(0);
  PdePartCutInterface::operator=(pdepartwithfemtraditional);
  return *this;
}
std::string PdePartMultiplier::getClassName() const
{
  return "PdePartMultiplier";
}
/*--------------------------------------------------------------------------*/
const solver_options::pdepart::opts PdePartMultiplier::setOptions()
{
  return solver_options::pdepart::none;
}
/*--------------------------------------------------------------------------*/
void PdePartMultiplier::additionCouplings(alat::Matrix<alat::SparsityPatternSoft>& sparsitypatternsoft)const
{
  alat::SparsityPatternSoft& sparsitypatternLU = sparsitypatternsoft(_ivars[1],_ivars[0]);
  alat::SparsityPatternSoft& sparsitypatternUL = sparsitypatternsoft(_ivars[0],_ivars[1]);

  for(int il=0; il < indicesofcutedge.n_cols; il++)
  {
    int i0 = indicesofcutedge(0,il);
    int i1 = indicesofcutedge(1,il);
    int i2 = indicesofcutedge(2,il);
    int i3 = indicesofcutedge(3,il);

    sparsitypatternLU[il].insert( i0 );
    sparsitypatternLU[il].insert( i1 );
    sparsitypatternLU[il].insert( i2 );
    sparsitypatternLU[il].insert( i3 );

    sparsitypatternUL[i0].insert( il );
    sparsitypatternUL[i1].insert( il );
    sparsitypatternUL[i2].insert( il );
    sparsitypatternUL[i3].insert( il );
  }
}
/*--------------------------------------------------------------------------*/
void PdePartMultiplier::computeMatrixGlobal(alat::MatrixAllVariables& A, const alat::VectorAllVariables& u)const
{
  alat::MatrixOneVariable* Aul = dynamic_cast<alat::MatrixOneVariable*>(A.get(_ivars[0],_ivars[1])); assert(Aul);
  alat::MatrixOneVariable* Alu = dynamic_cast<alat::MatrixOneVariable*>(A.get(_ivars[1],_ivars[0])); assert(Alu);
  const alat::SparsityPattern* Sul = Aul->getSparsityPattern();
  alat::armavec* Vul = Aul->getValues();
  const alat::SparsityPattern* Slu = Alu->getSparsityPattern();
  alat::armavec* Vlu = Alu->getValues();
  for(int il=0; il < indicesofcutedge.n_cols; il++)
  {
    for(int ii=0;ii<4;ii++)
    {
      int i = indicesofcutedge(ii,il);
      double d = coefsofcutedge(ii, il);
      int posstart = Slu->rowstart(il);
      int posend = Slu->rowstop(il);
      bool found=false;
      for(int pos = posstart; pos < posend; pos++)
      {
        if(Slu->col(pos)==i)
        {
          (*Vlu)[pos] += d;
          found=true;
          break;
        }
      }
      assert(found);
      posstart = Sul->rowstart(i);
      posend = Sul->rowstop(i);
      found=false;
      for(int pos = posstart; pos < posend; pos++)
      {
        if(Sul->col(pos)==il)
        {
          (*Vul)[pos] += d;
          found=true;
          break;
        }
      }
      assert(found);
    }
  }
}
/*--------------------------------------------------------------------------*/
void PdePartMultiplier::computeResidualGlobal(alat::VectorAllVariables& r, const alat::VectorAllVariables& u)const
{
  alat::VectorOneVariable* ru = dynamic_cast<alat::VectorOneVariable*>(r.get(_ivars[0])); assert(ru);
  alat::VectorOneVariable* rl = dynamic_cast<alat::VectorOneVariable*>(r.get(_ivars[1])); assert(rl);
  const alat::VectorOneVariable* uu = dynamic_cast<const alat::VectorOneVariable*>(u.get(_ivars[0])); assert(ru);
  const alat::VectorOneVariable* ul = dynamic_cast<const alat::VectorOneVariable*>(u.get(_ivars[1])); assert(rl);
  for(int il=0; il < indicesofcutedge.n_cols; il++)
  {
    int i0 = indicesofcutedge(0,il);
    int i1 = indicesofcutedge(1,il);
    int i2 = indicesofcutedge(2,il);
    int i3 = indicesofcutedge(3,il);

    double d0 = coefsofcutedge(0, il);
    double d1 = coefsofcutedge(1, il);
    double d2 = coefsofcutedge(2, il);
    double d3 = coefsofcutedge(3, il);

    (*rl)[il] -= d0*(*uu)[i0] + d1*(*uu)[i1] + d2*(*uu)[i2] + d3*(*uu)[i3];
    (*ru)[i0] -= d0*(*ul)[il];
    (*ru)[i1] -= d1*(*ul)[il];
    (*ru)[i2] -= d2*(*ul)[il];
    (*ru)[i3] -= d3*(*ul)[il];
  }
}

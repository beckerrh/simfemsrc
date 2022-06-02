#include  "Alat/vectorallvariables.hpp"
#include  "Alat/strings.hpp"
#include  "Alat/stringvector.hpp"
#include  <cassert>
#include  <cmath>
#include  <sstream>

using namespace alat;

/*--------------------------------------------------------------------------*/
VectorAllVariables::~VectorAllVariables() {}
VectorAllVariables::VectorAllVariables() : alat::Vector<std::shared_ptr<alat::VectorOneVariableInterface> >() {}
VectorAllVariables::VectorAllVariables(int nvars) : alat::Vector<std::shared_ptr<alat::VectorOneVariableInterface> >()
{
  set_size(nvars);
}
VectorAllVariables::VectorAllVariables( const VectorAllVariables& vectorallvariables) : alat::Vector<std::shared_ptr<alat::VectorOneVariableInterface> >(vectorallvariables)
{
  set_size(vectorallvariables.size());
  for(int i=0; i < size(); i++)
  {
    (*this)[i] = vectorallvariables.get(i)->clone();
  }
}
VectorAllVariables& VectorAllVariables::operator=( const VectorAllVariables& vectorallvariables)
{
  assert(0);
  return *this;
}
std::string VectorAllVariables::getClassName() const
{
  std::stringstream ss;
  ss << "VectorAllVariables_" << this->size();
  return ss.str();
}
VectorAllVariables* VectorAllVariables::clone() const
{
  return new VectorAllVariables(*this);
}
const alat::VectorOneVariableInterface* VectorAllVariables::get(int i) const
{
  return (*this)[i].get();
}
alat::VectorOneVariableInterface* VectorAllVariables::get(int i)
{
  return (*this)[i].get();
}

/*--------------------------------------------------------------------------*/
std::ostream& alat::operator<<(std::ostream& os, const VectorAllVariables& v)
{
  v.save(os, arma::arma_ascii);
  // os << v.getClassName() << "\n";
  // for(int i=0; i < v.size(); i++)
  // {
  //   v.get(i)->writeAscii(os);
  //   os << "\n";
  // }
  return os;
}

/*--------------------------------------------------------------------------*/
void VectorAllVariables::loadhdf5(const std::string& filename, const alat::StringVector& names)
{
  assert(names.size()==size());
  alat::armaivec ncomps, ns;
  ncomps.load(arma::hdf5_name(filename+".h5","/ncomps"));
  ns.load(arma::hdf5_name(filename+".h5","/ns"));
  for(int i=0;i<size();i++)
  {
    arma::hdf5_name spec(filename+".h5", names[i]);
    (*this)[i]->loadhdf5(spec);
    (*this)[i]->setNcompAndN(ncomps[i], ns[i]);
  }
}

/*--------------------------------------------------------------------------*/
void VectorAllVariables::savehdf5(const std::string& filename, const alat::StringVector& names) const
{
  assert(names.size()==size());
  assert(names.size()==size());
  alat::armaivec ncomps(size()), ns(size());
  for(int i=0;i<size();i++)
  {
    ncomps[i] = (*this)[i]->ncomp();
    ns[i] = (*this)[i]->n();
  }
  ncomps.save(arma::hdf5_name(filename+".h5","/ncomps"));
  ns.save(arma::hdf5_name(filename+".h5","/ns",arma::hdf5_opts::append));
  for(int i=0;i<size();i++)
  {
    // arma::hdf5_name spec(filename+".h5", names[i],arma::hdf5_opts::append+arma::hdf5_opts::trans);
    arma::hdf5_name spec(filename+".h5", names[i], arma::hdf5_opts::append);
    (*this)[i]->savehdf5(spec);
  }
}

/*--------------------------------------------------------------------------*/
void VectorAllVariables::save(std::ostream& os, arma::file_type datatype) const
{
  for(int i=0; i < this->size(); i++)
  {
    (*this)[i]->save(os, datatype);
  }
}
/*--------------------------------------------------------------------------*/
void VectorAllVariables::fillzeros()
{
  for(int i=0; i < this->size(); i++)
  {
    (*this)[i]->fillzeros();
  }
}
/*--------------------------------------------------------------------------*/
double VectorAllVariables::dot(const alat::VectorAllVariables* v) const
{
  double d = 0.0;
  for(int i=0; i < this->size(); i++)
  {
    d += (*this)[i]->dot( this->get(i));
  }
  return d;
}
double VectorAllVariables::norm() const {return sqrt( dot(this) );}
void VectorAllVariables::equal(double d)
{
  for(int i=0; i < this->size(); i++)
  {
    (*this)[i]->equal( d );
  }
}
void VectorAllVariables::equal(const alat::VectorAllVariables* v)
{
  for(int i=0; i < this->size(); i++)
  {
    (*this)[i]->equal( v->get(i) );
  }
}
void VectorAllVariables::add(const double& d, const alat::VectorAllVariables* v)
{
  for(int i=0; i < this->size(); i++)
  {
    (*this)[i]->add( d, v->get(i));
  }
}
/*--------------------------------------------------------------------------*/
// void VectorAllVariables::scalePerVariables(const alat::Vector<alat::armavec>& scales)
// {
//   for(int i=0; i < this->size(); i++)
//   {
//     (*this)[i]->scale(scales[i]);
//   }
// }
//
// void VectorAllVariables::scalePerVariablesInverse(const alat::Vector<alat::armavec>& scales)
// {
//   for(int i=0; i < this->size(); i++)
//   {
//    (*this)[i]->scaleinv(scales[i]);
//  }
// }

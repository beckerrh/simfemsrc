#include  "pdepart.hpp"

/*--------------------------------------------------------------------------*/
PdePart::~PdePart() {}
PdePart::PdePart(alat::StringList vars): solvers::PdePartWithIntegration(vars){}
PdePart::PdePart( const PdePart& pdepartwithfemtraditional): solvers::PdePartWithIntegration(pdepartwithfemtraditional)
{
  assert(0);
}
PdePart& PdePart::operator=( const PdePart& pdepartwithfemtraditional)
{
  assert(0);
  solvers::PdePartWithIntegration::operator=(pdepartwithfemtraditional);
  return *this;
}
std::string PdePart::getClassName() const
{
  return "PdePart";
}
double PdePart::_supg(double weight, double beta) const
{
  double h = pow(weight, 1.0/_meshinfo->dim);
  double a = beta;
  double b = 0.1*_localmodel->_diff/h;
  return _deltasupg*h/sqrt(a*a+b*b);
}

/*--------------------------------------------------------------------------*/
void PdePart::setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters)
{
  solvers::PdePartWithIntegration::setData(var2index, parameters);
  assert(_fems->size()==1);
  _deltasupg = parameters.doubles["deltasupg"];
  _gamma = parameters.doubles["gamma"];
  _symmetric = parameters.bools["symmetric"];
  _lumpedmass = parameters.bools["lumpedmass"];
  _ivar = 0;
  _ncomp = (*_fems)[_ivar]->getNcomp();
  _nlocal = (*_fems)[_ivar]->getNPerCell();
  _localmodel = dynamic_cast<const Model*>(_model);
  assert(_localmodel);
  // _beta.set_size(_mesh->getDimension());
  _beta.set_size(3);
  _beta.fill(0.0);
  const alat::VectorOneVariableInterface* datavector = _meshunit->getDataVector("beta");
  _betavec = dynamic_cast<const alat::VectorOneVariable*>(datavector);
  assert(_betavec);
  _femrt = dynamic_cast<const solvers::RT0*>(_meshunit->getFemData("beta"));
  assert(_femrt);
  assert(_localmodel);
  _betafct = _localmodel->getBeta().get();
  assert(_betafct);
}

#ifndef ___PdePart_hpp
#define ___PdePart_hpp

#include  "Solvers/pdepartwithintegration.hpp"
#include  "Solvers/rt0.hpp"
#include  "model.hpp"

/*--------------------------------------------------------------------------*/
class PdePart : public solvers::PdePartWithIntegration
{
protected:
  mutable alat::armavec _beta;
  const Model* _localmodel;
  int _ivar, _ncomp, _nlocal;
  double _deltasupg;
  double _supg(double weight, double beta) const;
  const alat::VectorOneVariable* _betavec;
  const solvers::RT0* _femrt;
  const solvers::InitialConditionInterface* _betafct;
  bool _lumpedmass;
  double _gamma;

public:
  ~PdePart();
  PdePart(alat::StringList vars);
  PdePart( const PdePart& pdepartwithfemtraditional);
  PdePart& operator=( const PdePart& pdepartwithfemtraditional);
  std::string getClassName() const;

  void setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters);
};

/*--------------------------------------------------------------------------*/
#endif

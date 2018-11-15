#include  "Perulangan/iterationinfo.hpp"
#include  <cassert>
#include  <cmath>

using namespace perulangan;

/*--------------------------------------------------------------------------*/
IterationInfoData::IterationInfoData() : maxiter(20), printstep(1), rtol(1e-6), gtol(1e-12), divfactor(1000.0) {}

/*--------------------------------------------------------------------------*/
IterationInfo::~IterationInfo(){}
IterationInfo::IterationInfo() : _firstresidual(0.0), _reductionrate(1.0), _iteration(0), _lastresidual(1.0), _printbuff("\n"), _printbufff("\n"){}
IterationInfo::IterationInfo( const IterationInfo& iterationinfo)
{
  _firstresidual = iterationinfo._firstresidual;
  _reductionrate = iterationinfo._reductionrate;
  _iteration = iterationinfo._iteration;
  _lastresidual = iterationinfo._lastresidual;
  _printbuff = iterationinfo._printbuff;
  _printbufff = iterationinfo._printbufff;
}
IterationInfo& IterationInfo::operator=( const IterationInfo& iterationinfo)
{
  assert(0);
  return *this;
}
std::string IterationInfo::getClassName() const
{
  return "IterationInfo";
}
IterationInfo* IterationInfo::clone() const
{
  return new IterationInfo(*this);
}

/*-------------------------------------------------------*/
std::ostream& IterationInfo::printLoopInformation(std::ostream& os) const
{
  os << "id/maxiter/printstep/rtol/gtol: " << _id<<"/"<<_iterationinfodata.maxiter<<"/"<<_iterationinfodata.printstep<<"/"<<_iterationinfodata.rtol<<"/"<<_iterationinfodata.gtol;
  return os;
}

/*-------------------------------------------------------*/
std::ostream& perulangan::operator<<(std::ostream& stream, const IterationInfo& iterationinfo)
{
  const IterationInfoData& data = iterationinfo.getIterationInfoData();
  stream << "rtol gtol maxiter _iteration " << data.rtol << " " << data.gtol << " " << data.maxiter << " " << iterationinfo.getNumberOfIterations();
  return stream;
}

/*-------------------------------------------------------*/
void IterationInfo::reset() const{_iteration = 0;}
const IterationInfoData& IterationInfo::getIterationInfoData() const{return _iterationinfodata;}
IterationInfoData& IterationInfo::getIterationInfoData() {return _iterationinfodata;}
double IterationInfo::getLastReductionRate() const{return _reductionrate;}
double IterationInfo::getPreviousReductionRate() const{return _previuousreductionrate;}
double IterationInfo::getMissingRate() const
{
  // ce qu'il reste a faire
  // soit le min entre gtol/res et rtol*firstres/res;
  return fmin(_iterationinfodata.gtol, _iterationinfodata.rtol*_firstresidual)/_lastresidual;
}
void IterationInfo::setId(std::string id){_id = id;}
const std::string& IterationInfo::getId() const {return _id;}
int IterationInfo::getNumberOfIterations() const{return _iteration;}
// const double& IterationInfo::getFirstResidual() const{return _firstresidual;}
const double& IterationInfo::getLastResidual() const{return _lastresidual;}
const double& IterationInfo::getGlobalTol() const {return _iterationinfodata.gtol;}
int IterationInfo::getMaxiter() const {return _iterationinfodata.maxiter;}

/*-------------------------------------------------------*/
void IterationInfo::printFirstIteration() const
{
  printf( "%3d %10s **=rate=** --- **%9.3e** (%6.1e %6.1e)%s", _iteration, getId().c_str(), _lastresidual, _firstresidual*_iterationinfodata.rtol, _iterationinfodata.gtol, _printbufff.c_str() );
  if(_firstresidual != _firstresidual)
  {
    assert(0);
  }
}
void IterationInfo::printIteration() const
{
  printf( "%3d %10s **%6.4f** --- **%9.3e**%s", _iteration, getId().c_str(), _reductionrate, _lastresidual, _printbuff.c_str() );
}

/*-------------------------------------------------------*/
void IterationInfo::checkIteration(perulanganEnums::iterationstatus& status, double res) const
{
  return checkIteration(status, res, _iteration, true);
}
void IterationInfo::checkIteration(perulanganEnums::iterationstatus& status, double residual, int iteration, bool print) const
{
  _iteration = iteration;
  if(iteration)
  {
    _previuousreductionrate = _reductionrate;
    _reductionrate = residual/_lastresidual;
  }
  // std::cerr << "IterationInfo::checkIteration() _iteration="<<_iteration<<" _lastresidual="<<_lastresidual<<" residual="<<residual << " print="<<print << "_iterationinfodata.printstep="<<_iterationinfodata.printstep<<" id"<<_id<<"\n";
  _lastresidual = residual;
  if(_iterationinfodata.printstep && !( iteration%_iterationinfodata.printstep ) && print && iteration)
  {
    printIteration();
  }
  if(iteration == 0)
  {
    status = perulanganEnums::IterationStatusRunning;
    // std::cerr << "_firstresidual="<<_firstresidual<< " " << getId() <<"\n";
    _firstresidual = residual;
    _lastresidual = residual;
    if(_iterationinfodata.printstep && !( iteration%_iterationinfodata.printstep ) && print)
    {
      printFirstIteration();
    }
  }
  // std::cerr << "IterationInfo::CheckIteration() " << getId() << " ? " << residual << " iteration=" << iteration<< " _iterationinfodata.maxiter=" << _iterationinfodata.maxiter << " status="<<perulanganEnums::iterationStatusToString(status)<<"\n";
  if( ( residual > _iterationinfodata.divfactor*_firstresidual )  ||  ( residual > _iterationinfodata.divfactor*_lastresidual ) )
  {
    if(print)
    {
      std::cerr << "IterationInfo::CheckIteration() BIG RESIDUAL _firstresidual _lastresidual residual "  << getId() << " ? iteration " << iteration << " _firstresidual " << _firstresidual << " _lastresidual " << _lastresidual << " residual " << residual << " _iterationinfodata.divfactor " << _iterationinfodata.divfactor << "\n";
      // assert(0);
    }
    status = perulanganEnums::IterationStatusDiverged;
    return;
  }
  else if(iteration >= _iterationinfodata.maxiter )
  {
    status = perulanganEnums::IterationStatusMaxIter;
    return;
  }
  _lastresidual = residual;
  bool ok = ( residual <= _iterationinfodata.rtol*_firstresidual )||( residual <= _iterationinfodata.gtol );
  if(ok)
  {
    status = perulanganEnums::IterationStatusConverged;
  }
  return;
}

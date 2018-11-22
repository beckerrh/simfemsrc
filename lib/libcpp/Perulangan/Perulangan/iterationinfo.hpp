#ifndef __Perulangan_IterationInfo_h
#define __Perulangan_IterationInfo_h

#include  <iostream>
#include  "enums.hpp"

/*--------------------------------------------------------------------------*/
namespace perulangan
{
  struct IterationInfoData
  {
    int maxiter, printstep;
    double rtol, gtol, divfactor;
    IterationInfoData();
  };
  class IterationInfo
  {
private:
    std::string _id;

    // user-defined parameters
    IterationInfoData _iterationinfodata;

    // information computed during execution
    mutable int _iteration;
    mutable double _lastresidual, _firstresidual, _reductionrate, _previuousreductionrate;

    void printFirstIteration() const;
    void printIteration() const;

public:
    ~IterationInfo();
    IterationInfo();
    IterationInfo( const IterationInfo& iterationinfo);
    IterationInfo& operator=( const IterationInfo& iterationinfo);
    std::string getClassName() const;
    IterationInfo* clone() const;

    void reset() const;
    std::ostream& printLoopInformation(std::ostream& os) const;

    void setId(std::string id);
    const IterationInfoData& getIterationInfoData() const;
    IterationInfoData& getIterationInfoData();
    const std::string& getId() const;
    int getNumberOfIterations() const;
    const double& getLastResidual() const;
    const double& getGlobalTol() const;
    int getMaxiter() const;
    void checkIteration(perulanganEnums::iterationstatus& status, double res) const;
    void checkIteration(perulanganEnums::iterationstatus& status, double res, int iter, bool print = 1) const;
    double getLastReductionRate() const;
    double getPreviousReductionRate() const;
    double getMissingRate() const;

    mutable std::string _printbuff, _printbufff;
  };
  std::ostream& operator<<(std::ostream& stream, const IterationInfo& A);
}

/*--------------------------------------------------------------------------*/

#endif

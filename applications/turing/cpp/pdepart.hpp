#ifndef ___PdePart_hpp
#define ___PdePart_hpp

#include  "Solvers/pdepartp1.hpp"

/*--------------------------------------------------------------------------*/
  class PdePart : public solvers::PdePartP1
  {
  protected:
    alat::armavec _k;
    double _a, _b;
    void brussel(arma::subview_col<double> f, const arma::subview_col<double> u) const;
    void brussel_d(arma::mat& df, const arma::subview_col<double> u)const;

  public:
    ~PdePart();
    PdePart(alat::StringList vars);
    PdePart( const PdePart& pdepart);
    PdePart& operator=( const PdePart& pdepart);
    std::string getClassName() const;

    void setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters);

    void computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeMatrixCell(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const;
  };

/*--------------------------------------------------------------------------*/
#endif

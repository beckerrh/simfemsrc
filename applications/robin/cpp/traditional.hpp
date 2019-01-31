#ifndef ___Traditional_hpp
#define ___Traditional_hpp

#include  "pdepart.hpp"
#include  "model.hpp"

/*--------------------------------------------------------------------------*/
  namespace alat
  {
    class Node;
  }
  class Traditional : public PdePart
  {
  protected:
    bool _lumpedmass;

  public:
    ~Traditional();
    Traditional(alat::StringList vars);
    Traditional( const Traditional& pdepart);
    Traditional& operator=( const Traditional& pdepart);
    std::string getClassName() const;

    void setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters);
    const solver_options::pdepart::opts setOptions();
    void computeRhsBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeResidualBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeMatrixBdry(int color, int iK, int iS, int iil, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const;
    void computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeMatrixCell(int iK, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const;
  };

/*--------------------------------------------------------------------------*/
#endif

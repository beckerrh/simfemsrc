#ifndef ___NewNitsche_hpp
#define ___NewNitsche_hpp

#include  "traditional.hpp"

/*--------------------------------------------------------------------------*/
  class NewNitsche : public Traditional
  {
  protected:
    mutable arma::uvec _isI;
    void _setisI(int iK)const;
    const arma::uvec* _dofisbdry;

  public:
    ~NewNitsche();
    NewNitsche(alat::StringList vars);
    NewNitsche( const NewNitsche& nitschenew);
    NewNitsche& operator=( const NewNitsche& nitschenew);
    std::string getClassName() const;

    void setData(const alat::Map<std::string, int>& var2index, const solvers::Parameters& parameters);
    const solver_options::pdepart::opts setOptions();
    void computeRhsCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeResidualCell(int iK, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeMatrixCell(int iK, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const;
    // void computeRhsBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    // void computeResidualBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    // void computeMatrixBdry(int color, int iK, int iS, int iil, solvers::PdePartData::mat& mat, solvers::PdePartData::imat& mat_i, solvers::PdePartData::imat& mat_j, const solvers::PdePartData::vec& uloc)const;
  };

/*--------------------------------------------------------------------------*/
#endif

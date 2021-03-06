#ifndef ___Nitsche_hpp
#define ___Nitsche_hpp

#include  "traditional.hpp"

/*--------------------------------------------------------------------------*/
  class Nitsche : public Traditional
  {
  protected:

  public:
    ~Nitsche();
    Nitsche(alat::StringList vars);
    Nitsche( const Nitsche& nitsche);
    Nitsche& operator=( const Nitsche& nitsche);
    std::string getClassName() const;

    const solver_options::pdepart::opts setOptions();
    void computeRhsBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeResidualBdry(int color, int iK, int iS, int iil, solvers::PdePartData::vec& floc, const solvers::PdePartData::vec& uloc)const;
    void computeMatrixBdry(int color, int iK, int iS, int iil, solvers::PdePartData::mat& mat, const solvers::PdePartData::vec& uloc)const;
  };

/*--------------------------------------------------------------------------*/
#endif

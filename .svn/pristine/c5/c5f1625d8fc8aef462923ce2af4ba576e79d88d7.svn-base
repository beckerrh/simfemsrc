#ifndef ___P1Cut_hpp
#define ___P1Cut_hpp

#include  "p1hansbo.hpp"

/*--------------------------------------------------------------------------*/
class P1Cut : public P1Hansbo
{
protected:

public:
  ~P1Cut();
  P1Cut();
  P1Cut( const P1Cut& P1cut);
  P1Cut& operator=( const P1Cut& P1cut);
  std::string getClassName() const;

  arma::vec beta;
  arma::mat Iin, Iex, Id, P, Q, Pnc, Qnc;

  void initFem(int ivar, const mesh::MeshUnitInterface* mesh, const solvers::MeshInfo* meshinfo, int ncomp);
  void setCell(int iK);
  void computeBeta(int iK);
};

/*--------------------------------------------------------------------------*/
#endif

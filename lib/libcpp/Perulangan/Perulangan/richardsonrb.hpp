#ifndef __Perulangan_RichardsonRB_h
#define __Perulangan_RichardsonRB_h

#include  "iterativesolverwithpreconditioner.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class RichardsonRB : public IterativeSolverWithPreconditioner
  {
protected:
    mutable double _rnorm;
    mutable int _niterafterrestar;
    int _nvectors, _nvars, _nshift;
    int _nextupdate, _nextproduct;
    mutable int _nmemory, _nextmemory;
    std::string _type, _solutiontype;
    mutable const alat::GhostMatrix* _ghostmatrix;
    mutable double _condition, _conditionmax, _conditionmean;
    mutable arma::mat _H;
    mutable alat::armavec _b, _x;
    alat::GhostVector& getV(int i) const;
    alat::GhostVector& getAV(int i, int ivar = 0) const;
    void _solveSmallSystem(perulanganEnums::iterationstatus& status) const;
    virtual void _computeSmallSystem(int index, int nmemory) const = 0;
    virtual void _addvector(alat::GhostVector& u, int nmemory) const = 0;
    virtual void _matrixVectorProduct(int index) const = 0;

public:
    ~RichardsonRB();
    RichardsonRB(const std::string& type, int nvectors, const std::string& solutiontype);
    RichardsonRB( const RichardsonRB& richardsonoptimal);
    RichardsonRB& operator=( const RichardsonRB& richardsonoptimal);
    std::string getClassName() const;

    void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
    void addUpdate(perulanganEnums::iterationstatus& status, const alat::GhostVector& w, alat::GhostVector& u) const;
    void restart();
  };
}

/*--------------------------------------------------------------------------*/

#endif

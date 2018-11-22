#ifndef __Alat_UmfMatrix_h
#define __Alat_UmfMatrix_h

#include  <string>
#include  "Alat/matrixonevariableinterface.hpp"

/*--------------------------------------------------------------------------*/

namespace alat
{
  class SparsityPattern;
}

namespace alat
{
  class UmfMatrix
  {
private:
    const alat::MatrixOneVariableInterface* _sparsematrix;

protected:
    double* Control;
    double* Info;
    void* Symbolic, * Numeric;

public:
    ~UmfMatrix();
    UmfMatrix();
    UmfMatrix( const UmfMatrix& umfmatrixbase);
    UmfMatrix& operator=( const UmfMatrix& umfmatrixbase);
    std::string getClassName() const;
    UmfMatrix* clone() const;

    void save(std::ostream& os, arma::file_type datatype = arma::arma_binary) const;
    void reInit(const alat::MatrixOneVariableInterface* sparsematrix);
    void computeLu();
    void solve(arma::vec& x, const arma::vec& b) const;
    void solveTranspose(arma::vec& x, const arma::vec& b) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

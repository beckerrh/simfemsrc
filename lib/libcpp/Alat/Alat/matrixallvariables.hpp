#ifndef __Alat_MatrixAllVariables_h
#define __Alat_MatrixAllVariables_h

#include  "Alat/matrix.hpp"
#include  "Alat/matrixonevariableinterface.hpp"
#include  "Alat/matrixonevariable.hpp"
#include  "Alat/umfmatrix.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class StringSet;
  class SystemAssembleMatrix;
  class VectorAllVariables;

  struct MatrixInOne
  {
    MatrixOneVariable matrix;
    alat::armaivec offsets;
    mutable alat::armavec out, in;
  };

  class MatrixAllVariables : public alat::Matrix<std::shared_ptr<alat::MatrixOneVariableInterface> >
  {
protected:
    // MatrixOneVariable _sparsematrix;
    // alat::UmfMatrix _umf;
    // alat::armaivec _offsets;
    // mutable alat::armavec _out, _in;

public:
    ~MatrixAllVariables();
    MatrixAllVariables();
    MatrixAllVariables(int nvars, int mvars);
    MatrixAllVariables( const MatrixAllVariables& matrixallvariables);
    MatrixAllVariables& operator=( const MatrixAllVariables& matrixallvariables);
    std::string getClassName() const;
    MatrixAllVariables* clone() const;

    const alat::MatrixOneVariableInterface* get(int i, int j) const;
    alat::MatrixOneVariableInterface* get(int i, int j);
    void fillzeros();
    void save(std::ostream& os, arma::file_type datatype = arma::arma_binary) const;
    void matrixVectorProduct(alat::VectorAllVariables* out, const alat::VectorAllVariables* in, double d = 1.0) const;
    void addMatrix(const MatrixAllVariables& matrix, double d = 1.0);
    // void solve(alat::VectorAllVariables& out, const alat::VectorAllVariables& in);
    void reInit(MatrixInOne& matrixinone) const;
    void compute(MatrixInOne& matrixinone) const;
  };
  std::ostream& operator<<(std::ostream& s, const MatrixAllVariables& matrix);
}

/*--------------------------------------------------------------------------*/

#endif

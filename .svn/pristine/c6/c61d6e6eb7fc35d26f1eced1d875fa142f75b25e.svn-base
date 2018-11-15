#ifndef __Alat_MatrixOneVariable_hpp
#define __Alat_MatrixOneVariable_hpp

#include  "Alat/matrixonevariableinterface.hpp"
#include  "Alat/sparsitypattern.hpp"
#include  "Alat/umfmatrix.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class MatrixOneVariable : public MatrixOneVariableInterface
  {
  protected:
    alat::SparsityPattern _sparsitypattern;
    arma::vec _values;
    alat::UmfMatrix _umf;

  public:
    ~MatrixOneVariable();
    MatrixOneVariable();
    MatrixOneVariable( const MatrixOneVariable& matrixonevariable);
    MatrixOneVariable& operator=( const MatrixOneVariable& matrixonevariable);
    std::string getClassName() const;
    MatrixOneVariable* clone() const;

    const alat::SparsityPattern* getSparsityPattern() const;
    const arma::vec* getValues() const;
    arma::vec* getValues();

    bool needsConnectivity() const;
    void set_size(int n, int m);
    void initSparsityPattern(const SparsityPatternSoft& sparsitypatternsoft);
    void save(std::ostream& os, arma::file_type datatype = arma::arma_binary) const;
    void write(std::ostream& os) const;
    void fillzeros();
    void matrixVectorProduct(alat::VectorOneVariableInterface* out, const alat::VectorOneVariableInterface* in, double d = 1.0) const;
    void addMatrix(const MatrixOneVariableInterface* matrix, double d = 1.0);
    void assemble(const arma::vec& Alocal, const alat::armaivec& indicesi, const alat::armaivec& indicesj);
    void rowIdentity(int index);
    void rowZero(int index);
    void solve(VectorOneVariableInterface* u, const VectorOneVariableInterface* f);
    void addEntriesForDirectSolver(int offsetivar, int offsetjvar, alat::SparsityPatternSoft& sparsitypatternsoft) const;
    void addMatrixForDirectSolver(int offsetivar, int offsetjvar, arma::vec& matrixvalues, const alat::SparsityPattern* sparsitypattern) const;
  };
  std::ostream& operator<<(std::ostream& os, const MatrixOneVariable& A);
}

/*--------------------------------------------------------------------------*/
#endif

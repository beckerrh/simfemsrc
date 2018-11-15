#ifndef __Alat_MatrixOneVariableArma_h
#define __Alat_MatrixOneVariableArma_h

#include  "Alat/armadillo.hpp"
#include  "Alat/matrixonevariableinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class MatrixOneVariableArma : public virtual arma::sp_mat, public virtual alat::MatrixOneVariableInterface
  {
  protected:
  public:
    ~MatrixOneVariableArma();
    MatrixOneVariableArma();
    MatrixOneVariableArma( const MatrixOneVariableArma& variablematrix);
    MatrixOneVariableArma& operator=( const MatrixOneVariableArma& variablematrix);
    std::string getClassName() const;
    MatrixOneVariableArma* clone() const;

    bool needsConnectivity() const;
    void fillzeros();
    void solve(VectorOneVariableInterface* u, const VectorOneVariableInterface* f);
    void set_size(int n, int m);
    void assemble(const arma::vec& Alocal, const alat::armaivec& indicesi, const alat::armaivec& indicesj);
    void rowIdentity(int index);
    void save(std::ostream& os, arma::file_type datatype = arma::arma_binary) const;
    void matrixVectorProduct(alat::VectorOneVariableInterface* out, const alat::VectorOneVariableInterface* in, double d = 1.0) const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

#ifndef __Alat_MatrixOneVariableInterface_hpp
#define __Alat_MatrixOneVariableInterface_hpp

#include  "Alat/interfacebase.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class SparsityPattern;
  class SparsityPatternSoft;
  class VectorOneVariableInterface;

  class MatrixOneVariableInterface : public virtual alat::InterfaceBase
  {
  public:
    ~MatrixOneVariableInterface();
    MatrixOneVariableInterface();
    MatrixOneVariableInterface( const MatrixOneVariableInterface& matrixonevariableinterface);
    MatrixOneVariableInterface& operator=( const MatrixOneVariableInterface& matrixonevariableinterface);
    std::string getClassName() const;
    MatrixOneVariableInterface* clone() const;

    virtual bool needsConnectivity() const=0;
    virtual void set_size(int n, int m)=0;
    virtual void save(std::ostream& os, arma::file_type datatype = arma::arma_binary) const=0;
    virtual void write(std::ostream& os) const;
    virtual void fillzeros()=0;
    virtual void matrixVectorProduct(alat::VectorOneVariableInterface* out, const alat::VectorOneVariableInterface* in, double d = 1.0) const;
    virtual void addMatrix(const MatrixOneVariableInterface* matrix, double d = 1.0);
    virtual void assemble(const alat::armamat& Alocal, const alat::armaivec& indicesi, const alat::armaivec& indicesj)=0;
    virtual void rowIdentity(int index);
    virtual void rowZero(int index);
    virtual void solve(VectorOneVariableInterface* u, const VectorOneVariableInterface* f);
    virtual void initSparsityPattern(const SparsityPatternSoft& sparsitypatternsoft);
    virtual const alat::SparsityPattern* getSparsityPattern() const;
    virtual const alat::armavec* getValues() const;
    virtual alat::armavec* getValues();
    virtual void addEntriesForDirectSolver(int offsetivar, int offsetjvar, alat::SparsityPatternSoft& sparsitypatternsoft) const;
    virtual void addMatrixForDirectSolver(int offsetivar, int offsetjvar, alat::armavec& matrixvalues, const alat::SparsityPattern* sparsitypattern) const;
 };
}

/*--------------------------------------------------------------------------*/
#endif

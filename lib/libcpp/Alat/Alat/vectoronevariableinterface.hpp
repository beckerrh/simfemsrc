#ifndef __Alat_VectorOneVariableInterface_hpp
#define __Alat_VectorOneVariableInterface_hpp

#include  "Alat/interfacebase.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class VectorOneVariable;

  class VectorOneVariableInterface : public virtual alat::InterfaceBase
  {
  private:
  protected:
  public:
    ~VectorOneVariableInterface();
    VectorOneVariableInterface();
    VectorOneVariableInterface( const VectorOneVariableInterface& vectoronevariableinterface);
    VectorOneVariableInterface& operator=( const VectorOneVariableInterface& vectoronevariableinterface);
    std::string getClassName() const;
    virtual std::unique_ptr<VectorOneVariableInterface> clone() const=0;

    virtual void loadhdf5(const arma::hdf5_name& spec)=0;
    virtual void savehdf5(const arma::hdf5_name& spec) const=0;
    virtual void save(std::ostream& os, arma::file_type datatype = arma::arma_binary) const=0;
    virtual int size() const=0;
    virtual void set_size(int n)=0;
    virtual void setNcompAndN(int ncomp, int n)=0;
    virtual int n() const=0;
    virtual int ncomp() const=0;
    virtual void scale(double s);
    virtual void fillzeros()=0;
    virtual double norm() const;
    virtual double dot(const alat::VectorOneVariableInterface* v) const=0;
    virtual void equal(const alat::VectorOneVariableInterface* v);
    virtual void equal(double d);
    virtual void add(const double& d, const alat::VectorOneVariableInterface* v);
    virtual void setVectorFromDirectSolver(int offset, const arma::vec& u);
    virtual void addVectorRhsForDirectSolver(int offset, arma::vec& f) const;
    virtual void scaleIntVector(const alat::armaivec& count)=0;

    virtual void assemble(const alat::armaimat& indices, const arma::mat& local, double d=1.0)=0;
    virtual void extract(const alat::armaimat& indices, arma::mat& local) const=0;
  };
}

/*--------------------------------------------------------------------------*/
#endif

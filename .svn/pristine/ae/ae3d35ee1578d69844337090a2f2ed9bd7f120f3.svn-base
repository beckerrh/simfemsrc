#ifndef __Alat_VectorAllVariables_h
#define __Alat_VectorAllVariables_h

#include  "Alat/vector.hpp"
#include  "Alat/vectoronevariableinterface.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  class StringSet;
  class StringVector;

  class VectorAllVariables : public alat::Vector<std::shared_ptr<alat::VectorOneVariableInterface> >
  {
public:
    ~VectorAllVariables();
    VectorAllVariables();
    VectorAllVariables(int nvars);
    VectorAllVariables( const VectorAllVariables& vectorallvariables);
    VectorAllVariables& operator=( const VectorAllVariables& vectorallvariables);
    std::string getClassName() const;
    VectorAllVariables* clone() const;

    const alat::VectorOneVariableInterface* get(int i) const;
    alat::VectorOneVariableInterface* get(int i);

    void fillzeros();
    double norm() const;
    double dot(const alat::VectorAllVariables* v) const;
    void equal(const alat::VectorAllVariables* v);
    void equal(double d);
    void add(const double& d, const alat::VectorAllVariables* v);
    // void scalePerVariables(const alat::Vector<arma::vec>& scales);
    // void scalePerVariablesInverse(const alat::Vector<arma::vec>& scales);
    void loadhdf5(const std::string& filename, const alat::StringVector& names);
    void savehdf5(const std::string& filename, const alat::StringVector& names) const;
    void save(std::ostream& os, arma::file_type datatype = arma::arma_binary) const;
   };
  std::ostream& operator<<(std::ostream& s, const VectorAllVariables& v);
}

/*--------------------------------------------------------------------------*/

#endif

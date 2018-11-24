#ifndef __Solvers_IntegrationFormulaInterface_h
#define __Solvers_IntegrationFormulaInterface_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/node.hpp"
#include  "Alat/vector.hpp"
#include  "Alat/armadillo.hpp"

/*--------------------------------------------------------------------------*/
namespace solvers
{
  class IntegrationFormulaInterface : public alat::InterfaceBase
  {
protected:
    alat::Vector<alat::Node> _c;
    alat::armavec _w;
    std::string getInterfaceName() const;

public:
    ~IntegrationFormulaInterface();
    IntegrationFormulaInterface(int n);
    IntegrationFormulaInterface( const IntegrationFormulaInterface& integrationformulainterface);
    IntegrationFormulaInterface& operator=( const IntegrationFormulaInterface& integrationformulainterface);
    std::string getClassName() const;

    virtual int n() const;
    virtual double weight(int k) const;
    virtual const alat::Node& point(int k) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

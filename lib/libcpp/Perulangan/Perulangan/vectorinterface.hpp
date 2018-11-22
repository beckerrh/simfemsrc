#ifndef __Perulangan_VectorInterface_h
#define __Perulangan_VectorInterface_h

#include  "Alat/interfacebase.hpp"

/*--------------------------------------------------------------------------*/
namespace perulangan
{
  class VectorInterface : public virtual alat::InterfaceBase
  {
protected:
    std::string getInterfaceName() const;

public:
    ~VectorInterface();
    VectorInterface();
    VectorInterface( const perulangan::VectorInterface& vectorinterface);
    VectorInterface& operator=( const perulangan::VectorInterface& vectorinterface);
    virtual VectorInterface* clone() const;

    virtual double norm() const;
    virtual double scalarProduct(const perulangan::VectorInterface* v) const;
    virtual void equal(const perulangan::VectorInterface* v);
    virtual void equal(double d);
    virtual void add(const double& d, const perulangan::VectorInterface* v);
    virtual void scale(const double& d);
    virtual void fillzeros();
    virtual std::ostream& writeAscii(std::ostream& os) const=0;
  };
}

/*--------------------------------------------------------------------------*/

#endif

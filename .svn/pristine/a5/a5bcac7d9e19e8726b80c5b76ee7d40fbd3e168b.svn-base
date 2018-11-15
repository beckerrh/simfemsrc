#ifndef __Alat_InterfaceBase_h
#define __Alat_InterfaceBase_h

#include  <string>

/*--------------------------------------------------------------------------*/

namespace alat
{
  template<class T>
  class Set;
}
namespace alat
{
  class InterfaceBase
  {
protected:
    mutable int _debug_level;
    virtual std::string getInterfaceName() const;
    void _notWritten(std::string function = "") const;
    void _message_string(std::string function, std::string message) const;
    void _warning_string(std::string function, std::string message) const;
    void _error_string(std::string function, std::string message) const;
    void _error_string(std::string function, std::string message, std::string value) const;
    void _error_string(std::string function, std::string message, int value) const;
    void _error_string(std::string function, std::string message, int value, int value1) const;
    void _error_string(std::string function, std::string message, const alat::Set<std::string>& values) const;

public:
    virtual ~InterfaceBase()=0;// {}
    InterfaceBase();
    InterfaceBase( const InterfaceBase& interfacebase);
    InterfaceBase& operator=( const InterfaceBase& interfacebase);
    void setDebugLevel(int debug_level) const;
    virtual std::string getClassName() const = 0;
  };
  inline InterfaceBase::~InterfaceBase() {}
}

/*--------------------------------------------------------------------------*/

#endif

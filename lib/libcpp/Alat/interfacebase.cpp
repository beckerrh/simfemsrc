#include  "Alat/interfacebase.hpp"
#include  "Alat/set.hpp"
#include  <cassert>
#include  <iostream>
#include  <stdlib.h>

using namespace alat;

/*--------------------------------------------------------------------------*/
InterfaceBase::InterfaceBase(): _debug_level(0)
{
    // std::cerr << "InterfaceBase::InterfaceBase()\n";
}
InterfaceBase::InterfaceBase( const InterfaceBase& interfacebase):_debug_level(interfacebase._debug_level){}
InterfaceBase& InterfaceBase::operator=( const InterfaceBase& interfacebase)
{
  _debug_level = interfacebase._debug_level;
  return *this;
}
std::string InterfaceBase::getInterfaceName() const
{
  return getClassName();
}
void InterfaceBase::setDebugLevel(int debug_level) const {_debug_level=debug_level;}

/*--------------------------------------------------------------------------*/
void InterfaceBase::_notWritten(std::string function) const
{
  std::cerr<<"*** ERROR in "<<getClassName()<<std::endl;
  std::cerr<<"\""<<getInterfaceName()<<"::"<<function<<"()\" not written !! ";
  std::cerr<<"(last class : \""<<getClassName()<<"\")\n";
  assert(0);
  exit(1);
}

/*--------------------------------------------------------------------------*/
void InterfaceBase::_message_string(std::string function, std::string message) const
{
  std::cerr<< "*** in "<< getClassName() << "::" << function<<"(): " << message  << "\n";
}

/*--------------------------------------------------------------------------*/
void InterfaceBase::_warning_string(std::string function, std::string message) const
{
  std::cerr<< "*** WARNING in "<< getClassName() << "::" << function<<"(): " << message  << "\n";
}

/*--------------------------------------------------------------------------*/
void InterfaceBase::_error_string(std::string function, std::string message) const
{
  std::cerr<< "*** ERROR in "<< getClassName() << "::" << function<<"(): " << message  << "\n";
  assert(0);
  exit(1);
}

/*--------------------------------------------------------------------------*/
void InterfaceBase::_error_string(std::string function, std::string message, std::string value) const
{
  std::cerr<< "*** ERROR in "<< getClassName() << "::" << function<<"(): "  << message << " \"" << value  << "\"\n";
  assert(0);
  exit(1);
}
void InterfaceBase::_error_string(std::string function, std::string message, int value) const
{
  std::cerr<< "*** ERROR in "<< getClassName() << "::" << function<<"(): " << message << " \"" << value  << "\"\n";
  assert(0);
  exit(1);
}
void InterfaceBase::_error_string(std::string function, std::string message, int value, int value2) const
{
  std::cerr<< "*** ERROR in "<< getClassName() << "::" << function<<"(): " << message << " \"" << value << " \" "<< " \"" << value2  << "\"\n";
  assert(0);
  exit(1);
}
void InterfaceBase::_error_string(std::string function, std::string message, const alat::Set<std::string>& values) const
{
  std::cerr<< "*** ERROR in "<< getClassName() << "::" << function<<"(): " << message << " " << values  << "\n";
  assert(0);
  exit(1);
}

#include  "Alat/ghost.hpp"
#include  <iostream>

using namespace alat;

/*--------------------------------------------------------------------------*/
Ghost::~Ghost(){}
Ghost::Ghost() : _name("none"), _type("none"){}
Ghost::Ghost(const std::string name) : _name(name), _type("none") {}
Ghost::Ghost(const std::string name, const std::string type) : _name(name), _type(type) {}
Ghost::Ghost( const Ghost& ghost) : _name(ghost._name), _type(ghost._type)
{}
Ghost& Ghost::operator=( const Ghost& ghost)
{
  _name = ghost._name;
  _type = ghost._type;
  return *this;
}

std::string Ghost::getClassName() const
{
  return "Ghost";
}

std::ostream& alat::operator<<(std::ostream& os, const Ghost& g)
{
  os << "(Name/type:) " << g.getName() <<"/"<< g.getType();
  return os;
}

/*--------------------------------------------------------------------------*/
void Ghost::setName(const std::string& name) {_name = name;}
const std::string& Ghost::getName() const{return _name;}

void Ghost::setType(const std::string& type){_type = type;}
const std::string& Ghost::getType() const{return _type;}

/*--------------------------------------------------------------------------*/
bool Ghost::operator<(const Ghost& v) const
{
  if( getName() < v.getName() )
  {
    return true;
  }
  if( getName() > v.getName() )
  {
    return false;
  }
  if( getType() < v.getType() )
  {
    return true;
  }
  if( getType() > v.getType() )
  {
    return false;
  }
  return false;
}

bool Ghost::operator==(const Ghost& v) const
{
  return ( getType() == v.getType() ) && ( getName() == v.getName() );
}

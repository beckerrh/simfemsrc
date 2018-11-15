#include  "Alat/ghostvector.hpp"
#include  <iostream>

using namespace alat;

/*--------------------------------------------------------------------------*/
GhostVector::~GhostVector(){}
GhostVector::GhostVector() : Ghost(){}
GhostVector::GhostVector(const std::string name) : Ghost(name) {}
GhostVector::GhostVector(const std::string name, const std::string type) : Ghost(name,type) {}
GhostVector::GhostVector( const GhostVector& ghost) : Ghost(ghost){}
GhostVector& GhostVector::operator=( const GhostVector& ghost)
{
  Ghost::operator=(ghost);
  return *this;
}

std::string GhostVector::getClassName() const
{
  return "GhostVector";
}

std::ostream& alat::operator<<(std::ostream& os, const GhostVector& g)
{
  os << "(Name/type:) " << g.getName() <<"/"<< g.getType();
  return os;
}

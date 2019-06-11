#include  "Alat/ghostmatrix.hpp"
#include  <iostream>

using namespace alat;

/*--------------------------------------------------------------------------*/
GhostMatrix::~GhostMatrix(){}
GhostMatrix::GhostMatrix() : Ghost(){}
GhostMatrix::GhostMatrix(const std::string name) : Ghost(name) {}
GhostMatrix::GhostMatrix(const std::string name, const std::string type) : Ghost(name,type) {}
GhostMatrix::GhostMatrix( const GhostMatrix& ghost) : Ghost(ghost){}
GhostMatrix& GhostMatrix::operator=( const GhostMatrix& ghost)
{
  Ghost::operator=(ghost);
  return *this;
}

std::string GhostMatrix::getClassName() const
{
  return "GhostMatrix";
}

std::ostream& alat::operator<<(std::ostream& os, const GhostMatrix& g)
{
  os << "(Name/type:) " << g.getClassName() <<"/"<< g.getType();
  return os;
}

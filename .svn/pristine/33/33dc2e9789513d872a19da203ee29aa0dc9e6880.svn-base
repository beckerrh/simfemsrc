#include  "Solvers/variable.hpp"
#include  "Solvers/p1.hpp"
#include  <cassert>
#include  <iostream>

using namespace solvers;

/*--------------------------------------------------------------------------*/
Variable::~Variable() {}
Variable::Variable() {}
Variable::Variable( const Variable& variable)
{
	(*this).operator=(variable);
}
Variable::Variable(std::string name, int ncomp, solverEnums::fem::femtype fem) : _name(name), _ncomp(ncomp), _fem(fem){}

Variable& Variable::operator=( const Variable& variable)
{
  _name = variable._name;
  _ncomp = variable._ncomp;
  _fem = variable._fem;
  return *this;
}
std::string Variable::getClassName() const
{
  return "Variable";
}
Variable* Variable::clone() const
{
  return new Variable(*this);
}

/*--------------------------------------------------------------------------*/
std::string Variable::getName() const
{
	return _name;
}
int Variable::getNcomp() const
{
	return _ncomp;
}
solverEnums::fem::femtype Variable::getFem() const
{
	return _fem;
}

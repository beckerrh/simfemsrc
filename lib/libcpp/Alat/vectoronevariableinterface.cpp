#include  "Alat/vectoronevariableinterface.hpp"
#include  <cassert>

using namespace alat;

/*--------------------------------------------------------------------------*/
VectorOneVariableInterface::~VectorOneVariableInterface() {}
VectorOneVariableInterface::VectorOneVariableInterface(): alat::InterfaceBase(){}
VectorOneVariableInterface::VectorOneVariableInterface( const VectorOneVariableInterface& vectoronevariableinterface): alat::InterfaceBase(vectoronevariableinterface) {}
VectorOneVariableInterface& VectorOneVariableInterface::operator=( const VectorOneVariableInterface& vectoronevariableinterface)
{
	assert(0);
	alat::InterfaceBase::operator=(vectoronevariableinterface);
	return *this;
}
std::string VectorOneVariableInterface::getClassName() const
{
	return "VectorOneVariableInterface";
}

/*--------------------------------------------------------------------------*/
void VectorOneVariableInterface::scale(double s){_notWritten("scale");}
// void VectorOneVariableInterface::scale(const alat::armavec& scale){_notWritten("scale");}
// void VectorOneVariableInterface::scaleinv(const alat::armavec& scale){_notWritten("scaleinv");}
double VectorOneVariableInterface::norm() const{_notWritten("norm");return 0.0;}
// double VectorOneVariableInterface::scalarProduct(const alat::VectorOneVariableInterface* v) const{_notWritten("scalarProduct");return 0.0;}
void VectorOneVariableInterface::equal(const alat::VectorOneVariableInterface* v){_notWritten("equal");}
void VectorOneVariableInterface::equal(double d){_notWritten("equal");}
void VectorOneVariableInterface::add(const double& d, const alat::VectorOneVariableInterface* v){_notWritten("add");}
void VectorOneVariableInterface::setVectorFromDirectSolver(int offset, const alat::armavec& u){_notWritten("setVectorFromDirectSolver");}
void VectorOneVariableInterface::addVectorRhsForDirectSolver(int offset, alat::armavec& f) const{_notWritten("addVectorRhsForDirectSolver");}

#include  "Alat/matrixonevariableinterface.hpp"
#include  <cassert>

using namespace alat;

/*--------------------------------------------------------------------------*/
MatrixOneVariableInterface::~MatrixOneVariableInterface() {}
MatrixOneVariableInterface::MatrixOneVariableInterface(): alat::InterfaceBase(){}
MatrixOneVariableInterface::MatrixOneVariableInterface( const MatrixOneVariableInterface& matrixonevariableinterface): alat::InterfaceBase(matrixonevariableinterface)
{
	(*this).operator=(matrixonevariableinterface);
}
MatrixOneVariableInterface& MatrixOneVariableInterface::operator=( const MatrixOneVariableInterface& matrixonevariableinterface)
{
	assert(0);
	alat::InterfaceBase::operator=(matrixonevariableinterface);
	return *this;
}
std::string MatrixOneVariableInterface::getClassName() const
{
	return "MatrixOneVariableInterface";
}
MatrixOneVariableInterface* MatrixOneVariableInterface::clone() const
{
	assert(0);
	return NULL;
	// return new MatrixOneVariableInterface(*this);
}

void MatrixOneVariableInterface::write(std::ostream& os) const {_notWritten("write");}

/*--------------------------------------------------------------------------*/
void MatrixOneVariableInterface::matrixVectorProduct(alat::VectorOneVariableInterface* out, const alat::VectorOneVariableInterface* in, double d) const {_notWritten("matrixVectorProduct");}
void MatrixOneVariableInterface::addMatrix(const MatrixOneVariableInterface* matrix, double d) {_notWritten("addMatrix");}
void MatrixOneVariableInterface::rowIdentity(int index){_notWritten("rowIdentity");}
void MatrixOneVariableInterface::rowZero(int index){_notWritten("rowZero");}
void MatrixOneVariableInterface::solve(VectorOneVariableInterface* u, const VectorOneVariableInterface* f){_notWritten("solve");}
void MatrixOneVariableInterface::initSparsityPattern(const SparsityPatternSoft& sparsitypatternsoft)
{_notWritten("initSparsityPattern");}
const alat::SparsityPattern* MatrixOneVariableInterface::getSparsityPattern() const {_notWritten("getSparsityPattern"); return NULL;}
const arma::vec* MatrixOneVariableInterface::getValues()const {_notWritten("getValues");return NULL;}
arma::vec* MatrixOneVariableInterface::getValues() {_notWritten("getValues");return NULL;}
void MatrixOneVariableInterface::addEntriesForDirectSolver(int offsetivar, int offsetjvar, alat::SparsityPatternSoft& sparsitypatternsoft) const{_notWritten("addEntriesForDirectSolver");}
void MatrixOneVariableInterface::addMatrixForDirectSolver(int offsetivar, int offsetjvar, arma::vec& matrixvalues, const alat::SparsityPattern* sparsitypattern) const{_notWritten("addMatrixForDirectSolver");}

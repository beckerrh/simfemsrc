#include  "Alat/sparsitypattern.hpp"
#include  "Alat/umfmatrix.hpp"
#include  "stdlib.h"
#include  <cassert>
#include  <fstream>
#include  <iostream>

#define UMFPACK_OK       0
#define UMFPACK_INFO    90
#define UMFPACK_CONTROL 20

#define UMFPACK_A       ( 0 )     /* Ax=b		*/
#define UMFPACK_At      ( 1 )     /* A'x=b	*/

using namespace alat;

// --------------------------------------------------------------------------*/
// #ifdef ARMA_64BIT_WORD
// /*--------------------------------------------------------------------------*/
// 
// extern "C" long umfpack_dl_symbolic
// (
//   long n,
//   long m,
//   const long Ap [],
//   const long Ai [],
//   const double Ax [],
//   void** Symbolic,
//   const double Control [UMFPACK_CONTROL],
//   double Info [UMFPACK_INFO]
// );
// extern "C" long umfpack_dl_numeric
// (
//   const long Ap [],
//   const long Ai [],
//   const double Ax [],
//   void* Symbolic,
//   void** Numeric,
//   const double Control [UMFPACK_CONTROL],
//   double Info [UMFPACK_INFO]
// );
// extern "C" long umfpack_dl_solve
// (
//   long sys,
//   const long Ap [],
//   const long Ai [],
//   const double Ax [],
//   double X [],
//   const double B [],
//   void* Numeric,
//   const double Control [UMFPACK_CONTROL],
//   double Info [UMFPACK_INFO]
// );
// extern "C" void umfpack_dl_free_symbolic(void** Symbolic);
// extern "C" void umfpack_dl_free_numeric(void** Numeric);
// extern "C" void umfpack_dl_report_status(const double Control [UMFPACK_CONTROL], long status);
// extern "C" void umfpack_dl_report_info(const double Control [UMFPACK_CONTROL], const double Info [UMFPACK_INFO]);
// extern "C" long umfpack_dl_report_numeric(const char name [], void* Numeric, const double Control [UMFPACK_CONTROL]);
// extern "C" void umfpack_dl_defaults(const double Control [UMFPACK_CONTROL]);
// 
// /*--------------------------------------------------------------------------*/
// UmfMatrix::~UmfMatrix()
// {
//   umfpack_dl_free_symbolic (&Symbolic);
//   umfpack_dl_free_numeric (&Numeric);
//   if(Control)
//   {
//     delete[] Control;
//     Control = NULL;
//   }
//   if(Info)
//   {
//     delete[] Info;
//     Info = NULL;
//   }
// }
// 
// UmfMatrix::UmfMatrix() : Control(NULL), Info(NULL), Symbolic(NULL), Numeric(NULL), _sparsematrix(NULL)
// {
//   std::cerr << "ARMA_64BIT_WORD\n";
//   assert(0);
//   Control = new double[UMFPACK_CONTROL];
//   umfpack_dl_defaults(Control);
//   Control[0] = 2;
// }
// 
// UmfMatrix::UmfMatrix( const UmfMatrix& umfmatrixbase) : Control(NULL), Info(NULL), Symbolic(NULL), Numeric(NULL)
// {
//   Control = new double[UMFPACK_CONTROL];
//   umfpack_dl_defaults(Control);
//   Control[0] = 2;
//   // assert(0);
// }
// 
// UmfMatrix& UmfMatrix::operator=( const UmfMatrix& umfmatrixbase)
// {
//   assert(0);
//   return *this;
// }
// 
// std::string UmfMatrix::getClassName() const
// {
//   return "UmfMatrix";
// }
// 
// UmfMatrix* UmfMatrix::clone() const
// {
//   assert(0);
//   return NULL;
// //return new UmfMatrix(*this);
// }
// 
// /*---------------------------------------------------------*/
// void UmfMatrix::save(std::ostream& os, arma::file_type datatype) const
// {
//   char name[] = "toto";
//   umfpack_dl_report_numeric(name, (void*) &Numeric, Control);
// }
// 
// /*-------------------------------------------------------------*/
// void UmfMatrix::reInit(const alat::MatrixOneVariableInterface* sparsematrix)
// {
//   _sparsematrix = sparsematrix;
//   umfpack_dl_free_symbolic (&Symbolic);
// 
//   const alat::SparsityPattern& sparsitypattern = *_sparsematrix->getSparsityPattern();
// 
//   long n = sparsitypattern.n();
//   const long* sb = &( *sparsitypattern.rowstart().begin() );
//   const long* cb = &( *sparsitypattern.col().begin() );
// 
//   long status = umfpack_dl_symbolic(n, n, sb, cb, NULL, &Symbolic, Control, Info);
// 
//   if(status != UMFPACK_OK)
//   {
//     umfpack_dl_report_info(Control, Info);
//     umfpack_dl_report_status(Control, status);
//     std::string datatype("ascii");
//     std::ofstream file("SPARSITYPATTERN_NOT_OK");
//     sparsitypattern.save(file, arma::arma_ascii);
//     std::cerr<<"*** ERROR UmfMatrix::reInit(: umfpack_symbolic failed\n";
//     assert(0);
//     exit(1);
//   }
// }
// 
// /*----------------------------------------------------------*/
// void UmfMatrix::computeLu()
// {
//   const alat::SparsityPattern& sparsitypattern = *_sparsematrix->getSparsityPattern();
//   const alat::armavec& mat = *_sparsematrix->getValues();
//   //
//   // construct LU
//   //
//   umfpack_dl_free_numeric (&Numeric);
//   const long* sb = &( *sparsitypattern.rowstart().begin() );
//   const long* cb = &( *sparsitypattern.col().begin() );
//   const double* mb = &( mat[0] );
//   long status = umfpack_dl_numeric(sb, cb, mb, Symbolic, &Numeric, Control, Info);
//   if(status != UMFPACK_OK)
//   {
//     umfpack_dl_report_info(Control, Info);
//     umfpack_dl_report_status(Control, status);
//     std::ofstream file("MATRIX_NOT_OK");
//     _sparsematrix->save(file, arma::arma_ascii);
//     _sparsematrix->save(std::cerr, arma::arma_ascii);
//     std::cerr<<"*** ERROR UmfMatrix::ComputeLu():umfpack_numeric failed\n";
//     assert(0);
//     exit(1);
//   }
//   //   umfpack_report_numeric("LU von A\n",Numeric,Control);
// }
// 
// /*----------------------------------------------------------*/
// void UmfMatrix::solve(alat::armavec& x, const alat::armavec& b) const
// {
//   const alat::SparsityPattern& sparsitypattern = *_sparsematrix->getSparsityPattern();
//   const alat::armavec& mat = *_sparsematrix->getValues();
//   assert( x.size() == b.size() );
//   assert( x.size() == sparsitypattern.n() );
//   const long* sb = &( *sparsitypattern.rowstart().begin() );
//   const long* cb = &( *sparsitypattern.col().begin() );
//   const double* mb = &( mat[0] );
//   double* xb = &( *x.begin() );
//   const double* bb = &( *b.begin() );
//   long status = umfpack_dl_solve (UMFPACK_At, sb, cb, mb, xb, bb, Numeric, Control, Info);
// 
//   if(status != UMFPACK_OK)
//   {
//     umfpack_dl_report_info(Control, Info);
//     umfpack_dl_report_status(Control, status);
//     std::ofstream file("MATRIX_NOT_OK");
//     _sparsematrix->save(file, arma::arma_ascii);
//     std::cerr<<"*** ERROR UmfMatrix::Solve(): umfpack_dl_solve failed\n";
//     assert(0);
//     exit(1);
//   }
// }
// 
// /*----------------------------------------------------------*/
// void UmfMatrix::solveTranspose(alat::armavec& x, const alat::armavec& b) const
// {
//   const alat::SparsityPattern& sparsitypattern = *_sparsematrix->getSparsityPattern();
//   const alat::armavec& mat = *_sparsematrix->getValues();
//   const long* sb = &( *sparsitypattern.rowstart().begin() );
//   const long* cb = &( *sparsitypattern.col().begin() );
//   const double* mb = &( mat[0] );
//   double* xb = &( *x.begin() );
//   const double* bb = &( *b.begin() );
//   long status = umfpack_dl_solve (UMFPACK_A, sb, cb, mb, xb, bb, Numeric, Control, Info);
// 
//   if(status != UMFPACK_OK)
//   {
//     umfpack_dl_report_info(Control, Info);
//     umfpack_dl_report_status(Control, status);
//     std::cerr<<"*** ERROR UmfMatrix::solveTranspose():umfpack_dl_solve failed\n";
//     assert(0);
//     exit(1);
//   }
// }
// 
// /*--------------------------------------------------------------------------*/
// #else
// /*--------------------------------------------------------------------------
extern "C" int umfpack_di_symbolic
(
  int n,
  int m,
  const int Ap [],
  const int Ai [],
  const double Ax [],
  void** Symbolic,
  const double Control [UMFPACK_CONTROL],
  double Info [UMFPACK_INFO]
);
extern "C" int umfpack_di_numeric
(
  const int Ap [],
  const int Ai [],
  const double Ax [],
  void* Symbolic,
  void** Numeric,
  const double Control [UMFPACK_CONTROL],
  double Info [UMFPACK_INFO]
);
extern "C" int umfpack_di_solve
(
  int sys,
  const int Ap [],
  const int Ai [],
  const double Ax [],
  double X [],
  const double B [],
  void* Numeric,
  const double Control [UMFPACK_CONTROL],
  double Info [UMFPACK_INFO]
);
extern "C" void umfpack_di_free_symbolic(void** Symbolic);
extern "C" void umfpack_di_free_numeric(void** Numeric);
extern "C" void umfpack_di_report_status(const double Control [UMFPACK_CONTROL], int status);
extern "C" void umfpack_di_report_info(const double Control [UMFPACK_CONTROL], const double Info [UMFPACK_INFO]);
extern "C" int umfpack_di_report_numeric(const char name [], void* Numeric, const double Control [UMFPACK_CONTROL]);
extern "C" void umfpack_di_defaults(const double Control [UMFPACK_CONTROL]);

/*--------------------------------------------------------------------------*/
UmfMatrix::~UmfMatrix()
{
  umfpack_di_free_symbolic (&Symbolic);
  umfpack_di_free_numeric (&Numeric);
  if(Control)
  {
    delete[] Control;
    Control = NULL;
  }
  if(Info)
  {
    delete[] Info;
    Info = NULL;
  }
}

UmfMatrix::UmfMatrix() : Control(NULL), Info(NULL), Symbolic(NULL), Numeric(NULL), _sparsematrix(NULL)
{
  Control = new double[UMFPACK_CONTROL];
  umfpack_di_defaults(Control);
  Control[0] = 2;
}

UmfMatrix::UmfMatrix( const UmfMatrix& umfmatrixbase) : Control(NULL), Info(NULL), Symbolic(NULL), Numeric(NULL)
{
  Control = new double[UMFPACK_CONTROL];
  umfpack_di_defaults(Control);
  Control[0] = 2;
}

UmfMatrix& UmfMatrix::operator=( const UmfMatrix& umfmatrixbase)
{
  assert(0);
  return *this;
}

std::string UmfMatrix::getClassName() const
{
  return "UmfMatrix";
}

UmfMatrix* UmfMatrix::clone() const
{
  assert(0);
  return NULL;
//return new UmfMatrix(*this);
}

/*---------------------------------------------------------*/
void UmfMatrix::save(std::ostream& os, arma::file_type datatype) const
{
  char name[] = "toto";
  umfpack_di_report_numeric(name, (void*) &Numeric, Control);
}

/*-------------------------------------------------------------*/
void UmfMatrix::reInit(const alat::MatrixOneVariableInterface* sparsematrix)
{
  _sparsematrix = sparsematrix;
  umfpack_di_free_symbolic (&Symbolic);

  const alat::SparsityPattern& sparsitypattern = *_sparsematrix->getSparsityPattern();

  int n = sparsitypattern.n();
  const int* sb = &( *sparsitypattern.rowstart().begin() );
  const int* cb = &( *sparsitypattern.col().begin() );

  int status = umfpack_di_symbolic(n, n, sb, cb, NULL, &Symbolic, Control, Info);

  if(status != UMFPACK_OK)
  {
    umfpack_di_report_info(Control, Info);
    umfpack_di_report_status(Control, status);
    std::string datatype("ascii");
    std::ofstream file("SPARSITYPATTERN_NOT_OK");
    sparsitypattern.save(file, arma::arma_ascii);
    // std::cerr<<"*** ERROR UmfMatrix::reInit(: umfpack_symbolic failed\n";
    assert(0);
    exit(1);
  }
}

/*----------------------------------------------------------*/
void UmfMatrix::computeLu()
{
  const alat::SparsityPattern& sparsitypattern = *_sparsematrix->getSparsityPattern();
  const alat::armavec& mat = *_sparsematrix->getValues();
  //
  // construct LU
  //
  umfpack_di_free_numeric (&Numeric);
  const int* sb = &( *sparsitypattern.rowstart().begin() );
  const int* cb = &( *sparsitypattern.col().begin() );
  const double* mb = &( mat[0] );
  int status = umfpack_di_numeric(sb, cb, mb, Symbolic, &Numeric, Control, Info);
  if(status != UMFPACK_OK)
  {
    umfpack_di_report_info(Control, Info);
    umfpack_di_report_status(Control, status);
    std::ofstream file("MATRIX_NOT_OK");
    _sparsematrix->write(file);
    // _sparsematrix->save(file, arma::arma_ascii);
    // _sparsematrix->save(std::cerr, arma::arma_ascii);
    // std::cerr<<"*** ERROR UmfMatrix::ComputeLu():umfpack_numeric failed\n";
    assert(0);
    exit(1);
  }
  //   umfpack_report_numeric("LU von A\n",Numeric,Control);
}

/*----------------------------------------------------------*/
void UmfMatrix::solve(alat::armavec& x, const alat::armavec& b) const
{
  const alat::SparsityPattern& sparsitypattern = *_sparsematrix->getSparsityPattern();
  const alat::armavec& mat = *_sparsematrix->getValues();
  assert( x.size() == b.size() );
  assert( x.size() == sparsitypattern.n() );
  const int* sb = &( *sparsitypattern.rowstart().begin() );
  const int* cb = &( *sparsitypattern.col().begin() );
  const double* mb = &( mat[0] );
  double* xb = &( *x.begin() );
  const double* bb = &( *b.begin() );
  int status = umfpack_di_solve (UMFPACK_At, sb, cb, mb, xb, bb, Numeric, Control, Info);

  if(status != UMFPACK_OK)
  {
    umfpack_di_report_info(Control, Info);
    umfpack_di_report_status(Control, status);
    std::ofstream file("MATRIX_NOT_OK");
    _sparsematrix->save(file, arma::arma_ascii);
    std::cerr<<"*** ERROR UmfMatrix::Solve(): umfpack_di_solve failed\n";
    assert(0);
    exit(1);
  }
}

/*----------------------------------------------------------*/
void UmfMatrix::solveTranspose(alat::armavec& x, const alat::armavec& b) const
{
  const alat::SparsityPattern& sparsitypattern = *_sparsematrix->getSparsityPattern();
  const alat::armavec& mat = *_sparsematrix->getValues();
  const int* sb = &( *sparsitypattern.rowstart().begin() );
  const int* cb = &( *sparsitypattern.col().begin() );
  const double* mb = &( mat[0] );
  double* xb = &( *x.begin() );
  const double* bb = &( *b.begin() );
  int status = umfpack_di_solve (UMFPACK_A, sb, cb, mb, xb, bb, Numeric, Control, Info);

  if(status != UMFPACK_OK)
  {
    umfpack_di_report_info(Control, Info);
    umfpack_di_report_status(Control, status);
    std::cerr<<"*** ERROR UmfMatrix::solveTranspose():umfpack_di_solve failed\n";
    assert(0);
    exit(1);
  }
}

/*--------------------------------------------------------------------------*/
//#endif
/*--------------------------------------------------------------------------*/

#undef UMFPACK_OK
#undef UMFPACK_INFO
#undef UMFPACK_CONTROL

#undef UMFPACK_A
#undef UMFPACK_At

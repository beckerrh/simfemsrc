#ifndef __Perulangan_TrivialPreconditioner_h
#define __Perulangan_TrivialPreconditioner_h

#include  "Perulangan/preconditioner.hpp"

/*--------------------------------------------------------------------------*/

namespace perulangan
{
  class TrivialPreconditioner : public perulangan::Preconditioner
  {
private:
protected:
public:
    ~TrivialPreconditioner();
    TrivialPreconditioner();
    TrivialPreconditioner( const TrivialPreconditioner& trivialpreconditioner);
    TrivialPreconditioner& operator=( const TrivialPreconditioner& trivialpreconditioner);
    std::string getClassName() const;
    TrivialPreconditioner* clone() const;

    int getNVectors() const;
    void reInit();
    void computePreconditioner();
    void solveApproximate(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f, int iteration) const;
    void setsmoothtype(std::string smoothtype);
  };
}

/*--------------------------------------------------------------------------*/

#endif

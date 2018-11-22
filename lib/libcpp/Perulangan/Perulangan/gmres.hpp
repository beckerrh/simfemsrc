#ifndef __Perulangan_Gmres_h
#define __Perulangan_Gmres_h

#include  "richardsonoptimal.hpp"

/*--------------------------------------------------------------------------*/
namespace perulangan
{
  class Gmres : public RichardsonOptimal
  {
public:
    ~Gmres();
    Gmres(const std::string& type, int nvectors, const std::string& solutiontype);
    Gmres( const Gmres& gmres);
    Gmres& operator=( const Gmres& gmres);
    std::string getClassName() const;

    int getNVectors() const;
    void solve(perulanganEnums::iterationstatus& status, const alat::GhostMatrix& A, alat::GhostVector& u, const alat::GhostVector& f) const;
  };
}

/*--------------------------------------------------------------------------*/

#endif

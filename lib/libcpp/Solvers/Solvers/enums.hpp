#ifndef  __Solvers_enum_h
#define  __Solvers_enum_h

#include  <string>

/*---------------------------------------------*/
namespace solverEnums
{
  enum assemblematrixtype {Diagonal, Full};
  namespace fem{
    enum femtype {None, Global, P1, P2, CR1, RT0, OWN};
    std::string femtypeToString(femtype f);
  }
}

#endif

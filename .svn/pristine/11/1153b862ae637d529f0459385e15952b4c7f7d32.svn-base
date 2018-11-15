#include  "Solvers/enums.hpp"
#include  <assert.h>
#include  <iostream>

namespace solverEnums
{
  namespace fem
  {
    std::string femtypeToString(femtype f)
    {
      if(f==None) {return "None";}
      else if(f==Global) {return "Global";}
      else if(f==P1) {return "P1";}
      else if(f==P2) {return "P2";}
      else if(f==CR1) {return "CR1";}
      else if(f==RT0) {return "RT0";}
      else if(f==OWN) {return "OWN";}
      std::cerr<<"unknown femtype " <<  f<<"\n"; assert(0); exit(1);
      return "";
    }
  }
}

#include  "simfem/timestepping.hpp"

// typedef double (*nlfct)(double u, double v);

//Hopf
// double global_Du = 0.0;
// double global_Dv = 0.0;
// double global_a = 1.0;
// double global_b = 2.1;
// double global_T = 50.0;
// double global_eps = 0.01;

//spatial pattern
// double global_Du = 0.001;
// double global_Dv = 0.01;
// double global_a = 1.0;
// double global_b = 2.0;
// double global_T = 30.0;
// double global_eps = 0.01;

//Turing
double global_Du = 0.0001;
double global_Dv = 0.01;
double global_a = 1.0;
double global_b = 1.5;
double global_T = 1.0;
double global_eps = 0.01;

/*---------------------------------------------------------------------------*/
int main(int argc, char** argv) 
{
  if (argc !=4)
  {
    printf("%s needs arguments <nx, nt, scheme>\n", argv[0]);
    exit(1);
  }    
  int nx = atoi(argv[1]);
  int nt = atoi(argv[2]);
  std::string scheme = argv[3];

  TimeSteppingData timesteppingdata;
  timesteppingdata.T = global_T;
  timesteppingdata.nt = nt;
  timesteppingdata.nx = nx;
  timesteppingdata.scheme = scheme;
  Nonlinearity nonlinearity("brusselator");
  nonlinearity.set_Du(global_Du);
  nonlinearity.set_Dv(global_Dv);
  nonlinearity.set_eps(global_eps);
  nonlinearity.set_a(global_a);
  nonlinearity.set_b(global_b);
  TimeStepping timestepping(timesteppingdata, nonlinearity);
  timestepping.run();
  return 0;
}
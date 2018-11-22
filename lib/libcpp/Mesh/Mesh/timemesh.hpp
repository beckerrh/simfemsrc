#ifndef  __Mesh_TimeMesh_h
#define  __Mesh_TimeMesh_h

#include  "Alat/interfacebase.hpp"
#include  "Alat/armadillo.hpp"

/*-----------------------------------------*/
namespace mesh
{
  struct TimeMeshData
  {
    int n;
    double first, last;
    TimeMeshData(int n_in=1, double first_in=0.0, double last_in=1.0) : n(n_in), first(first_in), last(last_in) {}
  };

  class TimeMesh : public virtual alat::InterfaceBase
  {
private:
    arma::vec _t;
    TimeMeshData _timemeshdata;

public:
    ~TimeMesh();
    TimeMesh();
    TimeMesh(const TimeMeshData timemeshdata);
    TimeMesh(const TimeMesh& timemesh);
    TimeMesh& operator=(const TimeMesh& timemesh);
    std::string getClassName() const;
    std::string getInfo() const;

    void setData(const TimeMeshData timemeshdata);
    int n() const;
    double t(int i) const;
    const arma::vec& t() const;
  };
}

#endif

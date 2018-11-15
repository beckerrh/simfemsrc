#ifndef  __MeshWrap_TimeMeshWrapper_h
#define  __MeshWrap_TimeMeshWrapper_h

#include  "Mesh/timemesh.hpp"

/*---------------------------------------------------------------------------*/
class TimeMeshWrapper : public mesh::TimeMesh
{
public:
    TimeMeshWrapper();
    TimeMeshWrapper(int n, double first, double last);
    void setData(int n, double first, double last);
};
std::ostream &operator<<(std::ostream &os, const TimeMeshWrapper& timemesh);

/*---------------------------------------------------------------------------*/
void wrapTimeMesh();


#endif

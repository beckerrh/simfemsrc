#include  "Mesh/timemesh.hpp"
#include  <fstream>
#include  <cassert>

using namespace mesh;

/*---------------------------------------------------------*/
TimeMesh::~TimeMesh() {}
TimeMesh::TimeMesh() : alat::InterfaceBase()
{
  setData(TimeMeshData());
}
TimeMesh::TimeMesh(const TimeMeshData timemeshdata) : alat::InterfaceBase()
{
	setData(timemeshdata);
}
TimeMesh::TimeMesh(const TimeMesh& timemesh) : alat::InterfaceBase(timemesh) {}
TimeMesh& TimeMesh::operator=(const TimeMesh& timemesh)
{
	assert(0);
	return *this;
}
std::string TimeMesh::getClassName() const
{
	return "mesh::TimeMesh";
}
std::string TimeMesh::getInfo() const
{
	std::stringstream ss;
	int n = _t.size();
	if(n==0)
	{
		_error_string("getInfo","n=",n);
	}
	ss << "first_last_n " << _t[0] << "_" << _t[n-1] << "_" << n-1;
	return ss.str();
}
void TimeMesh::setData(const TimeMeshData timemeshdata)
{
  _timemeshdata = timemeshdata;
  int n=_timemeshdata.n;
	_t.resize(n+1);
	double k = (_timemeshdata.last-_timemeshdata.first)/(double) n;
	_t[0] = _timemeshdata.first;
	for(int i=0; i<n; i++)
	{
		_t[i+1] = _t[i] + k;
	}
}

/*---------------------------------------------------------*/
int TimeMesh::n() const {return _t.size();}
double TimeMesh::t(int i) const{return _t[i];}
const alat::armavec& TimeMesh::t() const{return _t;}

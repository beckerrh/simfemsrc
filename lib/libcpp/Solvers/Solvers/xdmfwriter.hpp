#ifndef  __Solvers_XdmfWriter_h
#define  __Solvers_XdmfWriter_h

#include  "Mesh/timemesh.hpp"
#include  "Solvers/meshunitwithdatainterface.hpp"
#include  "Solvers/options.hpp"
#include  <string>

/*-----------------------------------------*/
namespace alat
{
  template<class T>
	class Vector;
}

namespace solvers
{
	class Variable;

  class XdmfWriter
  {
protected:
  mutable std::string _blank, _meshnametopo, _meshnamegeo;
  mutable std::stringstream _topo, _topodata, _geom, _geomdata;
  mutable int _ntimes, _nnodes;
  const solver_options::output_manager* _output_manager;
  const mesh::TimeMesh* _timemesh;
  const mesh::MeshUnitInterface* _mesh;
	std::string _meshfilenameh5, _solutionfilenameh5;

  void writeBegin(std::ofstream& file, int ntmax) const;
  void writeEnd(std::ofstream& file) const;
  void writeMesh(std::ofstream& file, int it) const;
  void writeScalars(std::ofstream& file, const std::string& solutionfilenameh5, int it, const alat::Vector<Variable>& vars) const;

public:
  XdmfWriter(const solver_options::output_manager* output_manager, const mesh::TimeMesh* timemesh, const mesh::MeshUnitInterface* mesh);
  void writeSolution(std::string filename, std::string name, const alat::Vector<Variable>& vars, int ntmax=-1) const;
  void writeData(std::string filename, const alat::Vector<Variable>& vars, int ntmax=-1) const;
  };
}

#endif

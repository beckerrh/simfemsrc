#include  "Alat/vector.hpp"
#include  "Mesh/meshunitinterface.hpp"
#include  "Solvers/xdmfwriter.hpp"
#include  "Solvers/variable.hpp"
#include  "Solvers/feminterface.hpp"
#include  <fstream>
#include  <iomanip>
#include  <sstream>

using namespace solvers;

/*---------------------------------------------------------*/
XdmfWriter::XdmfWriter(const solver_options::output_manager* output_manager, const mesh::TimeMesh* timemesh, const mesh::MeshUnitInterface* mesh) : _output_manager(output_manager), _timemesh(timemesh), _mesh(mesh)
{
  _blank = "  ";
}

/*---------------------------------------------------------*/
void XdmfWriter::writeSolution(std::string filename, std::string name, const alat::Vector<Variable>& vars, int ntmax) const
{
  std::ofstream file(filename.c_str());
  writeBegin(file, ntmax);
  for(int it=0;it<_ntimes;it++)
  {
    writeMesh(file, it);

    std::string solutionfilenameh5 = _output_manager->getSolutionFileName(meshEnums::Plain, name, it) + ".h5";
    writeScalars(file, solutionfilenameh5, it, vars);
  }
  writeEnd(file);
}
/*---------------------------------------------------------*/
void XdmfWriter::writeData(std::string filename, const alat::Vector<Variable>& vars, int ntmax) const
{
  std::ofstream file(filename.c_str());
  writeBegin(file, ntmax);
  for(int it=0;it<_ntimes;it++)
  {
    writeMesh(file, it);

    std::string solutionfilenameh5 = _output_manager->getDataFileName(meshEnums::Plain, it) + ".h5";
    writeScalars(file, solutionfilenameh5, it, vars);
  }
  writeEnd(file);
}

/*---------------------------------------------------------*/
void XdmfWriter::writeScalars(std::ofstream& file, const std::string& solutionfilenameh5, int it, const alat::Vector<Variable>& vars) const
{
  for(int ivar=0; ivar<vars.size(); ivar++)
  {
    std::string name = vars[ivar].getClassName();
    int ncomp =  vars[ivar].getNcomp();

    for(int icomp=0;icomp<ncomp;icomp++)
    {
      std::stringstream atr, atrdata, atrdatabase;
      atr << "<Attribute Name=\""<<name<< icomp << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
      atrdata << "<DataItem ItemType=\"HyperSlab\" Dimensions=\"1 " << _nnodes << "\" >\n";
      atrdata << "<DataItem Dimensions=\"3 2\"  Format=\"XML\"> \n";
      atrdata << 0 << " " << icomp*_nnodes << "\n" << 1 << " " << 1 << "\n" << 1 << " "<< _nnodes << "\n";
      atrdata << "</DataItem>\n";
      atrdata << "<DataItem Dimensions=\"1 " << ncomp*_nnodes;
      atrdata << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
      atrdatabase << solutionfilenameh5 << ":/"<<name<<"\n";

      file << std::setw (8) << _blank << atr.str();
      file << std::setw (10) << _blank << atrdata.str();
      file << std::setw (12) << _blank << atrdatabase.str();
      file << std::setw (10) << _blank << "</DataItem>\n";
      file << std::setw (10) << _blank << "</DataItem>\n";
      file << std::setw (8) << _blank << "</Attribute>\n";
    }
  }
  file << std::setw (6) << _blank << "</Grid>\n";
}
/*---------------------------------------------------------*/
void XdmfWriter::writeMesh(std::ofstream& file, int it) const
{
  if(it==0)
  {
    file << std::setw (6) << _blank << "<Grid Name=\"mesh\" GridType=\"Uniform\">\n";

    file << std::setw (8) << _blank << _topo.str();
    file << std::setw (10) << _blank << _topodata.str();
    file << std::setw (12) << _blank << _meshnametopo;
    file << std::setw (10) << _blank << "</DataItem>\n";
    file << std::setw (8) << _blank << "</Topology>\n";

    file << std::setw (8) << _blank << _geom.str();
    file << std::setw (10) << _blank << _geomdata.str();
    file << std::setw (12) << _blank << _meshnamegeo;
    file << std::setw (10) << _blank << "</DataItem>\n";
    file << std::setw (8) << _blank << "</Geometry>\n";
  }
  else
  {
    file << std::setw (6) << _blank << "<Grid>\n";
    std::string bizarre = "<xi:include xpointer=\"xpointer(//Grid[@Name=&quot;TimeSeries&quot;]/Grid[1]/*[self::Topology or self::Geometry])\" />";
    file << std::setw (8) << _blank << bizarre << "\n";
  }
  std::stringstream timevalue;
  timevalue << "<Time Value=\" " << _timemesh->t(it) << "\" />\n";
  file << std::setw (8) << _blank << timevalue.str();
}

/*---------------------------------------------------------*/
void XdmfWriter::writeBegin(std::ofstream& file, int ntmax) const
{
  int ncells = _mesh->getNCells();
  int nodesperelement = _mesh->getNNodesPerCell();
  _nnodes = _mesh->getNNodes();
  int dimension = 3;//_mesh->getDimension();
  std::string toptype = _mesh->getXdmfTopoType();

  std::string meshfilenameh5 = _output_manager->getMeshFileName() + ".h5";
  _meshnametopo = meshfilenameh5 + ":/PlainMesh/NodesAndNodesOfCells/nodes_of_cells\n";
  _meshnamegeo = meshfilenameh5 + ":/PlainMesh/NodesAndNodesOfCells/nodes\n";

  _topo << "<Topology NumberOfElements=\"" << ncells << "\" TopologyType=\"" << toptype <<"\" ";
  _topo << "NodesPerElement=\"" << nodesperelement <<"\">\n";
  _topodata << "<DataItem Dimensions=\"" << ncells << " " << nodesperelement;
  _topodata << "\" NumberType=\"UInt\" Format=\"HDF\">\n";
  _geom << "<Geometry GeometryType=\"XYZ\">\n";
  _geomdata << "<DataItem Dimensions=\"" << _nnodes  << " " << dimension;
  _geomdata << "\"  NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";

  std::string gridtime = "<Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";
  if(ntmax==-1)
  {
    _ntimes = _timemesh->n();
  }
  else
  {
    _ntimes = ntmax+1;
  }
  std::stringstream timeinfo1, timeinfo2, timeinfo3;
  timeinfo1 << "<Time TimeType=\"List\">" << "\n";
  timeinfo2 << "<DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"" << _ntimes << "\">" << "\n";
  timeinfo3 << std::setw(6);
  for(int it=0;it<_ntimes;it++) timeinfo3 << _timemesh->t(it)<< " ";
  timeinfo3 << "\n";
  file << "<?xml version=\"1.0\"?>\n";
  file << "<Xdmf Version=\"2.0\" xmlns:xi=\"http://www.w3.org/2001/XInclude\">\n";
  file << std::setw (2) << _blank << "<Domain>\n";
  file << std::setw (4) << _blank << gridtime;
}

/*---------------------------------------------------------*/
void XdmfWriter::writeEnd(std::ofstream& file) const
{
  file << std::setw (4) << _blank << "</Grid>\n";
  file << std::setw (2) << _blank << "</Domain>\n";
  file << std::setw (0) << "</Xdmf>\n";
  file.close();
}

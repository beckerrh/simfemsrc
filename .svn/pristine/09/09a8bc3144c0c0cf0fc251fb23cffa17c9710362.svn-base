#include  "Solvers/options.hpp"
#include  "Alat/directoryandfiles.hpp"
#include  <iomanip>
#include  <sstream>

using namespace solver_options;

/*--------------------------------------------------------------------------*/
std::string solver_options::output_manager_data_to_string(output_manager_data data)
{
  if(data==datadir) {return "datadir";}
  else if(data==solutionfilename) {return "solutionfilename";}
  else if(data==datafilename) {return "datafilename";}
  else if(data==meshfilename) {return "meshfilename";}
  else {assert(0); return "";}
}
std::ostream& solver_options::operator<<(std::ostream& os, const output_manager& om)
{
  for(alat::Map<output_manager_data, std::string>::const_iterator p=om.data.begin();p!=om.data.end();p++)
  {
    os << output_manager_data_to_string(p->first) << ": " << p->second << "\n";
  }
  return os;
}
void output_manager::init()
{
  path =   alat::_getPath() + "/";
}
//   // std::cerr << "_path="<<_path<<"\n";
//   std::string command;
//   if(alat::_directoryExists(data[datadir] + ".old"))
//   {
//     command = "rm -rf " + data[datadir] + ".old";
//     system( command.c_str() );
//   }
//   if(alat::_directoryExists(data[datadir]))
//   {
//     command = "mv -f "+data[datadir] + " " + data[datadir] + ".old";
//     system( command.c_str() );
//   }
//   command = "mkdir -p "+data[datadir];
//   system( command.c_str() );
//   command = "mkdir -p "+data[datadir] + "/Plain";
//   system( command.c_str() );
//   command = "mkdir -p "+data[datadir] + "/Boundary";
//   system( command.c_str() );
//   command = "mkdir -p "+data[datadir] + "/Interface";
//   system( command.c_str() );
// }
std::string output_manager::getSolutionFileName(meshEnums::meshunittype type, std::string name, int it) const
{
  if(it==-1) return path+data[datadir]+"/"+data[solutionfilename];
  std::stringstream ss;
  ss << "_" << std::setfill('0') << std::setw(6) << it;
  return path+data[datadir]+"/"+meshEnums::meshunitTypeToString(type)+"/"+data[solutionfilename]+name+ss.str();
}
std::string output_manager::getDataFileName(meshEnums::meshunittype type, int it) const
{
  if(it==-1) return path+data[datadir]+"/"+data[datafilename];
  std::stringstream ss;
  ss << "_" << std::setfill('0') << std::setw(6) << it;
  return path+data[datadir]+"/"+meshEnums::meshunitTypeToString(type)+"/"+data[datafilename]+ss.str();
}

std::string output_manager::getMeshFileName() const {return path+data[datadir]+"/"+data[meshfilename];}

/*--------------------------------------------------------------------------*/
bool solver_options::has_armamat(const opts& op) {return bool(op.flags & flag_armamat);}
bool solver_options::has_dynamic(const opts& op) {return bool(op.flags & flag_dynamic);}
bool solver_options::has_erasedatadir(const opts& op) {return bool(op.flags & flag_erasedatadir);}

std::ostream& solver_options::operator<<(std::ostream& os, const opts& op)
{
  os << "armamat(" << has_armamat(op)<< ") dynamic(" << has_dynamic(op)<< ") erasedatadir(" << has_erasedatadir(op)<<")";
  return os;
}
/*--------------------------------------------------------------------------*/
bool solver_options::errors::has_L2(const opts& op) {return bool(op.flags & flag_L2);}
bool solver_options::errors::has_H1(const opts& op) {return bool(op.flags & flag_H1);}
bool solver_options::errors::has_L1(const opts& op) {return bool(op.flags & flag_L1);}
bool solver_options::errors::has_Linf(const opts& op) {return bool(op.flags & flag_Linf);}
bool solver_options::errors::has_E(const opts& op) {return bool(op.flags & flag_E);}
std::ostream& solver_options::errors::operator<<(std::ostream& os, const solver_options::errors::opts& op)
{
  os << "L2(" << has_L2(op)<< ") Linf(" << has_Linf(op)<< ") L1(" << has_L1(op)<< ") H1(" << has_H1(op)<< ") E(" << has_E(op)<<")";
  return os;
}
/*--------------------------------------------------------------------------*/
bool solver_options::pdepart::has_cell(const opts& op) {return bool(op.flags & flag_cell);}
bool solver_options::pdepart::has_bdry(const opts& op) {return bool(op.flags & flag_bdry);}
bool solver_options::pdepart::has_iside(const opts& op) {return bool(op.flags & flag_iside);}
std::ostream& solver_options::pdepart::operator<<(std::ostream& os, const solver_options::pdepart::opts& op)
{
  os << "cell(" << has_cell(op)<< ") iside(" << has_iside(op)<< ") bdry(" << has_bdry(op)<< ")";
  return os;
}

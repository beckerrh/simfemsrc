#include  "Alat/chronometer.hpp"
#include  <cassert>
#include  <fstream>
#include  <iomanip>
#include  <iostream>

using namespace alat;

/*--------------------------------------------------------------------------*/
Chronometer::~Chronometer() {}
Chronometer::Chronometer() : _classname("") {}
Chronometer::Chronometer(const Chronometer& C) : _classname("") {}
Chronometer::Chronometer(std::string classname) : _classname(classname) {}
Chronometer& Chronometer::operator= (const Chronometer& C)
{
  assert(0);
  return *this;
}

/*--------------------------------------------------------------------------*/
void Chronometer::print(std::ostream& os) const
{
  double totaltime = total();
  os << "   Total        :  "  << getClassName() << "   " << std::setiosflags(std::ios::fixed)<< std::setprecision(4) << totaltime << " s\n";
  for(Chronometer::const_iterator p = _M.begin(); p != _M.end(); p++)
  {
    double singletime = p->second.read();
    if(_sum[p->first])
    {
      os << std::setiosflags(std::ios::left);
      os << std::setw(30) << p->first << "\t" << std::setiosflags(std::ios::fixed) << std::setprecision(5);
      os << std::resetiosflags(std::ios::left);
      os << std::setw(15) << singletime  << std::setw(12) << std::setiosflags(std::ios::fixed) << std::setprecision(2)<<  100.0*singletime/totaltime << "\% \n";
    }
    else
    {
      os << std::setiosflags(std::ios::left);
      os << std::setw(30) << p->first << "\t" << std::setiosflags(std::ios::fixed) << std::setprecision(5);
      os << std::resetiosflags(std::ios::left);
      os << std::setw(15) << singletime  << std::setw(12) << "---\n";
    }
  }
  os << "\n";
}

/*--------------------------------------------------------------------------*/
std::string Chronometer::getClassName() const
{
  return _classname;
}
void Chronometer::setClassName(const std::string& classname)
{
  _classname = classname;
}

/*--------------------------------------------------------------------------*/

const std::map<std::string, alat::StopWatch>& Chronometer::get() const
{
  return _M;
}
/*--------------------------------------------------------------------------*/
void Chronometer::enrol(std::string name, bool sum)
{
  _M[name];
  _sum[name] = sum;
}
/*--------------------------------------------------------------------------*/
void Chronometer::reset()
{
  for(Chronometer::const_iterator p = _M.begin(); p != _M.end(); p++)
  {
    reset(p->first);
  }
}

/*--------------------------------------------------------------------------*/
void Chronometer::reset(std::string name)
{
  get(name).reset();
}

/*--------------------------------------------------------------------------*/
void Chronometer::start(std::string name)
{
  if( get(name).running() )
  {
    std::cerr << "*** ERROR in Chronometer::start() already running function '" << name << "' in class "<< _classname <<"\n";
    assert(0);
    exit(1);
  }
  else
  {
    get(name).start();
  }
}

/*--------------------------------------------------------------------------*/
double Chronometer::stop(std::string name)
{
  return get(name).stop();
}

/*--------------------------------------------------------------------------*/
alat::StopWatch& Chronometer::get(std::string name)
{
  return _M[name];
}
/*--------------------------------------------------------------------------*/
double Chronometer::total() const
{
  double d = 0.;
  for(const_iterator p = _M.begin(); p != _M.end(); p++)
  {
    if( _sum[p->first] )
    {
      d += p->second.read();
    }
  }
  return d;
}

#ifndef __Alat_Chronometer_h
#define __Alat_Chronometer_h

#include  "Alat/map.hpp"
#include  "Alat/set.hpp"
#include  "stopwatch.hpp"
#include  <string>

/*--------------------------------------------------------------------------*/

namespace alat
{
  class Chronometer
  {
public:
    typedef alat::Map<std::string, alat::StopWatch>::const_iterator const_iterator;
    typedef alat::Map<std::string, alat::StopWatch>::iterator iterator;

private:
    std::string _classname;
    alat::Map<std::string, alat::StopWatch> _M;
    alat::Map<std::string, bool> _sum;

public:
    ~Chronometer();
    Chronometer();
    Chronometer(const Chronometer& C);
    Chronometer(std::string classname);
    Chronometer& operator= (const Chronometer& C);

    std::string getClassName() const;
    void setClassName(const std::string& classname);
    alat::StopWatch& get(std::string name);

    const std::map<std::string, alat::StopWatch>& get() const;
    void enrol(std::string name, bool sum = true);
    void reset();
    void reset(std::string name);
    void start(std::string name);
    double stop(std::string name);
    double total() const;
    void print(std::ostream& os) const;
  };
}

#endif

#ifndef __Alat_StopWatch_H
#define __Alat_StopWatch_H

#include  <ctime>

/*--------------------------------------------------------------------------*/

namespace alat
{
  class StopWatch
  {
private:
    bool _running;
    double _last_time;
    double _total;

    double seconds(void);

public:
    ~StopWatch();
    StopWatch();
    StopWatch(const StopWatch& s);
    StopWatch& operator=( const StopWatch& s );

    bool running();
    void reset();
    void start();
    double stop();
    double read() const;
  };
}

#endif

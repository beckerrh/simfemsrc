#ifndef __Alat_directoryandfiles_h
#define __Alat_directoryandfiles_h

#include  <string>

/*--------------------------------------------------------------------------*/

namespace alat
{
  std::string _getPath();
  bool _directoryExists(std::string strPath);
  bool _FileExists(std::string strPath);
}

/*--------------------------------------------------------------------------*/

#endif

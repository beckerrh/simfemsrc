#include  "Alat/directoryandfiles.hpp"
#include  <fstream>
#include  <iostream>
#include  <sys/stat.h>
#include  <sys/types.h>
#include  <dirent.h>
#include  <limits.h>
#include  <unistd.h>

namespace alat
{
  std::string _getPath()
  {
    char buffer[1024];
    return std::string(getcwd(buffer, sizeof(buffer)));
  }

  bool _directoryExists( std::string pzPath )
  {
    DIR* pDir;
    bool bExists = false;

    pDir = opendir ( pzPath.c_str() );

    if(pDir != NULL)
    {
      bExists = true;
      (void) closedir (pDir);
    }

    return bExists;
  }

/*--------------------------------------------------------------------------*/

  bool _FileExists(std::string strPath)
  {
    std::ifstream file( strPath.c_str() );
    bool filexists= file.is_open();
    file.close();
    return filexists;
  }
}

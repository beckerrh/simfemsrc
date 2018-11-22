#include  "Alat/getlinesplit.hpp"
#include  <iterator>
#include  <sstream>

using namespace alat;

/*---------------------------------------------------------*/

alat::StringVector alat::getLineSplit(std::ifstream& file)
{
  std::string line;
  while(line == "")
  {
    getline(file, line);
  }
  std::istringstream is(line);
  return alat::StringVector( std::istream_iterator<std::string>(is), std::istream_iterator<std::string>() );
}

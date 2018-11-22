#ifndef  __FadalightMesh_BoundaryInfo_h
#define  __FadalightMesh_BoundaryInfo_h

#include  "Alat/map.hpp"
#include  "Alat/armadillo.hpp"
#include  "Alat/vector.hpp"

/*---------------------------------------------------------*/
namespace FadalightMesh
{
  class BoundaryInfo
  {
protected:
    alat::IntMap  _index_of_col;
    alat::armaivec _colors;
    alat::Vector<alat::armaivec> _cells;
    alat::Vector<alat::armaivec> _sides;
    alat::Vector<alat::armaivec> _sideids_of_cells;

    void _ConstructIntOfColor();
    int _IndexOfColor(int color) const;

public: ~BoundaryInfo();
    BoundaryInfo();
    BoundaryInfo(const BoundaryInfo& boundaryinfo);
    BoundaryInfo(const std::string& filename);
    BoundaryInfo& operator=(const BoundaryInfo& boundaryinfo);

    void set_size(const std::map<int, int>& size_of_color);

    void read(const std::string& filename);
    void write(const std::string& filename, arma::file_type datatype = arma::arma_binary) const;

    // int getNBoundaries();
    int getNSides() const;
    alat::armaivec& getColors();
    const alat::armaivec& getColors() const;
    alat::armaivec& getCellsOfColor(int color);
    const alat::armaivec& getCellsOfColor(int color) const;
    alat::armaivec& getSidesOfColor(int color);
    const alat::armaivec& getSidesOfColor(int color) const;
    alat::armaivec& getSidesIdOfCellsOfColor(int color);
    const alat::armaivec& getSidesIdOfCellsOfColor(int color) const;
    // void removeColor(int color);
  };
}

#endif

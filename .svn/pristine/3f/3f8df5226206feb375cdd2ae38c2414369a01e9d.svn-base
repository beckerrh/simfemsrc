#ifndef __Alat_CubicSpline_hpp
#define __Alat_CubicSpline_hpp

#include  "armadillo.hpp"
#include  "node.hpp"

/*--------------------------------------------------------------------------*/
namespace alat
{
  struct CubicSplineNewtonInfo
  {
    arma::vec3 r, dr, ddr;
    double g, h, d;
  };
  std::ostream& operator<<(std::ostream& os, const CubicSplineNewtonInfo& info);

  class CubicSpline
  {
  protected:
    // armacube _C;
    arma::mat _A;
    arma::field<arma::vec3> _C;
    void bezier(arma::vec3& r, int i, double t) const;
    void dbezier(arma::vec3& dr, int i, double t) const;
    void ddbezier(arma::vec3& ddr, int i, double t) const;
    void _testsmooth() const;
    void _computeNewtonInfo(CubicSplineNewtonInfo& info, const arma::vec3& p, int i, double t) const;
    double _newtonMinDist(CubicSplineNewtonInfo& info, const arma::vec3& p, int i) const;
    double _newtonInterval(double t0, double t1, CubicSplineNewtonInfo& info, const arma::vec3& p, int i) const;

  public:
    ~CubicSpline();
    CubicSpline();
    CubicSpline( const CubicSpline& cubicspline);
    CubicSpline& operator=( const CubicSpline& cubicspline);
    std::string getClassName() const;

    void generate(const armamat& points);
    double signeddistance(const arma::vec3& p) const;
    void writeVtk(std::string filename) const;
  };
}

/*--------------------------------------------------------------------------*/
#endif

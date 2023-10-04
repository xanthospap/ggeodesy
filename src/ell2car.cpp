#include "geodesy.hpp"

void ell2car(double lambda, double phi, double h, const dso::Ellipsoid &e,
             double &x, double &y, double &z) noexcept {
  /* Eccentricity squared. */
  double e2{e.eccentricity_squared()};

  /* Radius of curvature in the prime vertical. */
  const double N{e.N(phi)};

  /* transform */
  core::geodetic2cartesian(lambda, phi, h, e2, N, x, y, z);

  /* Finished. */
  return;
}

Eigen::Matrix<double, 3, 1> dso::ell2car(const Eigen::Matrix<double, 3, 1> &lfh,
                                         const dso::Ellipsoid &e) noexcept {
  Eigen::Matrix<double, 3, 1> xyz;
  ell2car(lfh(0), lfh(1), lfh(2), e, xyz(0), xyz(1), xyz(2));
  return xyz;
}

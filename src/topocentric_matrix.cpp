#include "geodesy.hpp"
#include <cassert>

Eigen::Matrix<double, 3, 3> dso::topocentric_matrix(double lon,
                                                    double lat) noexcept {
  const double cf = std::cos(lat);
  const double sf = std::sin(lat);
  const double cl = std::cos(lon);
  const double sl = std::sin(lon);

  /* unit vector along east */
  const double e0 = -sl;
  const double e1 = cl;
  const double e2 = 0e0;

  /* unit vector along north */
  const double n0 = -sf * cl;
  const double n1 = -sf * sl;
  const double n2 = cf;

  /* unit vector along up */
  const double u0 = cf * cl;
  const double u1 = cf * sl;
  const double u2 = sf;

  /* Rotation matrix is: R = [n^t, e^T, u^T] */
  return (Eigen::Matrix<double, 3, 3>() << n0, n1, n2, e0, e1, e2, u0, u1, u2)
      .finished();
}

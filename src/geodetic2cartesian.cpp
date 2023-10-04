/* @file ell2car.cpp
 * @brief Transform geodetic (ellipsoidal) to cartesian coordinates.
 */

#include "geodesy.hpp"
#ifdef INLCUDE_GEODESY_CHECKS
#include <cassert>
#endif

void dso::ell2car(double lambda, double phi, double h, const dso::Ellipsoid &e,
                  double &x, double &y, double &z) noexcept {
#ifdef INCLUDE_GEODESY_CHECKS
  assert(dso::ellipsoidal_ok(lambda, phi));
#endif

  /* Eccentricity squared. */
  const double e2 = e.eccentricity_squared();

  /* Radius of curvature in the prime vertical. */
  const double N = e.N(phi);

  /* Trigonometric numbers. */
  const double sf = std::sin(phi);
  const double cf = std::cos(phi);
  const double sl = std::sin(lambda);
  const double cl = std::cos(lambda);

  /* Compute geocentric rectangular coordinates. */
  x = (N + h) * cf * cl;
  y = (N + h) * cf * sl;
  z = ((1e0 - e2) * N + h) * sf;

  /* Finished. */
  return;
}

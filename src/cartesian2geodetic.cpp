/* @file car2ell.cpp
 * @brief Transformation of ellipsoidal to cartesian coordinates.
 */

#include "geodesy.hpp"

/* @brief Cartesian to ellipsoidal/geodetic.
 *
 * Transform cartesian, geocentric coordinates (x, y, z) to ellipsoidal (i.e.
 * longitude, geodetic latitude and ellispoidal height), using the algorithm 
 * described in Fukushima, T., "Transformation from Cartesian to geodetic 
 * coordinates accelerated by Halley's method", J. Geodesy (2006)
 *
 * @param[in]  x      Cartesian x-component [m]
 * @param[in]  y      Cartesian y-component [m]
 * @param[in]  z      Cartesian z-component [m]
 * @param[in]  semi_major   Semi-major axis of the (reference) ellipsoid
 *                    (meters)
 * @param[in]  flattening   Flattening of the (reference) ellipsoid
 * @param[out] lambda Ellipsoidal longtitude in range (-π, π) [rad]
 * @param[out] phi    Ellipsoidal/geodetic latitude in range (-π/2, π/2) [rad]
 * @param[out] h      Ellipsoidal height [m]
 */
void dso::core::car2ell(double x, double y, double z, double semi_major,
                        double flattening, double &lambda, double &phi,
                        double &h) noexcept {
  /* Functions of ellipsoid parameters. */
  const double aeps2 = semi_major * semi_major * 1e-32;
  const double e2 = (2e0 - flattening) * flattening;
  const double e4t = e2 * e2 * 1.5e0;
  const double ep2 = 1e0 - e2;
  const double ep = std::sqrt(ep2);
  const double aep = semi_major * ep;

  /* Compute Coefficients of (Modified) Quartic Equation
   * Remark: Coefficients are rescaled by dividing by 'a'
   */

  /* Compute distance from polar axis squared. */
  const double p2=x * x + y * y;

  /* Compute longitude lambda. */
  if (p2) {
    lambda = std::atan2(y, x);
  } else {
    lambda = 0e0;
  }

  /* Ensure that Z-coordinate is unsigned. */
  const double absz = std::abs(z);

  /* Continue unless at the poles */
  if (p2 > aeps2) {
    /* Compute distance from polar axis. */
    const double p = std::sqrt(p2);
    /* Normalize. */
    const double s0 = absz / semi_major;
    const double pn = p / semi_major;
    const double zp = ep * s0;
    /* Prepare Newton correction factors. */
    const double c0 = ep * pn;
    const double c02 = c0 * c0;
    const double c03 = c02 * c0;
    const double s02 = s0 * s0;
    const double s03 = s02 * s0;
    const double a02 = c02 + s02;
    const double a0 = std::sqrt(a02);
    const double a03 = a02 * a0;
    const double d0 = zp * a03 + e2 * s03;
    const double f0 = pn * a03 - e2 * c03;
    /* Prepare Halley correction factor. */
    const double b0 = e4t * s02 * c02 * pn * (a0 - ep);
    const double s1 = d0 * f0 - b0 * s0;
    const double cp = ep * (f0 * f0 - b0 * c0);
    /* Evaluate latitude and height. */
    phi = std::atan(s1 / cp);
    const double s12 = s1 * s1;
    const double cp2 = cp * cp;
    h = (p * cp + absz * s1 - semi_major * std::sqrt(ep2 * s12 + cp2)) /
        std::sqrt(s12 + cp2);
  } else {
    /* Special case: pole. */
    phi = dso::DPI / 2e0;
    h = absz - aep;
  }

  /* Restore sign of latitude. */
  if (z < 0.e0)
    phi = -phi;

  return;
}

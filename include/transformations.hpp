/** @file
 * A list of frequently used geodetic functions.
 */

#ifndef __DSO_COORDINATE_TRANSFORMATIONS_HPP__
#define __DSO_COORDINATE_TRANSFORMATIONS_HPP__

#include "core/crd_transformations.hpp"
#include "core/crdtype_warppers.hpp"
#include "ellipsoid.hpp"

namespace dso {

/* @brief Geodetic (ellipsoidal) to cartesian coordinates.
 *
 * @tparam      E    The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @param[in]   lat  Geodetic latitude (-π/2, π/2) [rad]
 * @param[in]   lon  Geodetic longtitude in rage (-π, π) [rad]
 * @param[in]   h    Ellipsoidal height [m]
 * @param[out]  x    Cartesian x-component [m]
 * @param[out]  y    Cartesian y-component [m]
 * @param[out]  z    Cartesian z-component [m]
 */
template <ellipsoid E>
void geodetic2cartesian(double lat, double lon, double h, double &x, double &y,
                        double &z) noexcept {
  /* Eccentricity squared. */
  constexpr const double e2 = dso::eccentricity_squared<E>();

  /* Radius of curvature in the prime vertical. */
  const double Rn = dso::N<E>(lat);

  /* Trigonometric numbers. */
  const double sf = std::sin(lat);
  const double cf = std::cos(lat);
  const double sl = std::sin(lon);
  const double cl = std::cos(lon);

  /* Compute geocentric rectangular coordinates. */
  x = (Rn + h) * cf * cl;
  y = (Rn + h) * cf * sl;
  z = ((1e0 - e2) * Rn + h) * sf;

  /* Finished. */
  return;
}

/** @brief Cartesian to geodetic/ellipsoidal.
 *
 * Transform cartesian, geocentric coordinates (x, y, z) to ellipsoidal (i.e.
 * latitude, longtitude, ellispoidal height). All units are meters and
 * radians.
 * Fukushima, T., "Transformation from Cartesian to geodetic coordinates
 * accelerated by Halley's method", J. Geodesy (2006), 79(12): 689-693
 *
 * @tparam     E    The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @param[in]  x    Cartesian x-component (meters)
 * @param[in]  y    Cartesian y-component (meters)
 * @param[out] lat  Geodetic latitude (radians)
 * @param[out] lon  Geodetic longtitude (radians)
 * @param[out] h    Ellipsoidal height (meters)
 */
template <ellipsoid E>
void cartesian2geodetic(double x, double y, double z, double &lat, double &lon,
                        double &hgt) noexcept {
  /* Functions of ellipsoid parameters. */
  constexpr double a = ellipsoid_traits<E>::a;
  constexpr double aeps2 = a * a * 1e-32;
  constexpr double e2 = eccentricity_squared<E>();
  constexpr double e4t = e2 * e2 * 1.5e0;
  constexpr double ep2 = 1.0e0 - e2;
#if defined(__clang__)
  /* clang goes pedantic on this! std::sqrt is not constexpr */
  const double ep = std::sqrt(ep2);
  const double aep = a * ep;
#else
  constexpr double ep = std::sqrt(ep2);
  constexpr double aep = a * ep;
#endif

  /* Compute distance from polar axis squared. */
  const double p2 = x * x + y * y;

  /* Compute longitude. */
  lon = (p2) ? std::atan2(y, x) : 0e0;

  /* Ensure that Z-coordinate is unsigned. */
  const double absz = std::abs(z);

  if (p2 > aeps2) {
    /* Continue unless at the poles */
    /* Compute distance from polar axis. */
    const double p = std::sqrt(p2);
    /* Normalize. */
    const double s0 = absz / a;
    const double pn = p / a;
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
    lat = std::atan(s1 / cp);
    const double s12 = s1 * s1;
    const double cp2 = cp * cp;
    hgt = (p * cp + absz * s1 - a * std::sqrt(ep2 * s12 + cp2)) /
          std::sqrt(s12 + cp2);
  } else {
    /* Special case: pole. */
    lat = dso::DPI / 2e0;
    hgt = absz - aep;
  }

  /* Restore sign of latitude. */
  if (z < 0.e0)
    lat = -lat;

  /* Finished. */
  return;
}

template <typename C = CartesianCrd>
inline SphericalCrd cartesian2spherical(const /*CartesianCrd*/ C &v) noexcept {
  static_assert(dso::CoordinateTypeTraits<C>::isCartesian);
  SphericalCrd s;
  cartesian2spherical(v.x(), v.y(), v.z(), s.r(), s.lat(), s.lon());
  return s;
}

template <typename S = SphericalCrd>
inline CartesianCrd spherical2cartesian(const /*SphericalCrd*/ S &v) noexcept {
  static_assert(dso::CoordinateTypeTraits<S>::isSpherical);
  CartesianCrd s;
  spherical2cartesian(v.r(), v.lat(), v.lon(), s.x(), s.y(), s.z());
  return s;
}

template <ellipsoid E, typename G = GeodeticCrd>
CartesianCrd geodetic2cartesian(const /*GeodeticCrd*/ G &v) noexcept {
  static_assert(dso::CoordinateTypeTraits<G>::isGeodetic);
  CartesianCrd s;
  geodetic2cartesian<E>(v.lat(), v.lon(), v.hgt(), s.x(), s.y(), s.z());
  return s;
}

template <ellipsoid E, typename C = CartesianCrd>
GeodeticCrd cartesian2geodetic(const /*CartesianCrd*/ C &v) noexcept {
  static_assert(dso::CoordinateTypeTraits<C>::isCartesian);
  GeodeticCrd s;
  cartesian2geodetic<E>(v.x(), v.y(), v.z(), s.lat(), s.lon(), s.hgt());
  return s;
}
} /* namespace dso */

#endif

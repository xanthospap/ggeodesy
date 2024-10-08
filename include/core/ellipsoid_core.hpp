/** @file
 *  Core functions of ellipsoidal geometry. In geodesy we mostly use the
 *  semi-major axis a and the flattening f as defining paramters, hence here,
 *  we use these as fundamental parameters to compute derived geometric
 *  quantities.
 */

#ifndef __DSO_ELLIPSOID_GEOMETRY_CORE_HPP__
#define __DSO_ELLIPSOID_GEOMETRY_CORE_HPP__

#include <cmath>

namespace dso {

/** @brief core namespace holds the core of ellipsoid-related functions */
namespace core {

/** @brief Compute the squared eccentricity.
 *
 * Compute the squared (first) eccentricity (i.e. \f$e^2\f$) given the
 * flattening of an ellipsoid, aka \f$ e^2 = \frac{a^2-b^2}{a^2} = (2-f)*f \f$
 *
 * @param[in] f Flattening
 * @return Squared eccentricity
 */
inline constexpr double eccentricity_squared(double f) noexcept {
  return (2e0 - f) * f;
}

/** @brief Compute the third flattening
 *
 * Compute the third flattening (usually denoted as \f$ n \f$), given the
 * flattening \f$ f \f$, aka \f$ n = \frac{a-b}{a+b} = f/(2-f) \f$.
 *
 * @param[in] f Flattening
 * @return Third flattening \f$ n \f$
 */
inline constexpr double third_flattening(double f) noexcept {
  return f / (2e0 - f);
}

/** @brief Compute the semi-minor axis (aka \f$b\f$)
 *
 * Compute the semi-minor axis of an ellipsoid \f$b\f$), given the
 * flattening and the semi-major axis: \f$ \beta = \alpha * (1-f) \f$.
 * Units of semi-minor axis are the same as the input units of semi-major
 * axis (parameter a).
 *
 * @param[in] f Flattening
 * @param[in] a Semi-major axis
 * @return Semi-minor axis with units same as a.
 */
inline constexpr double semi_minor(double a, double f) noexcept {
  return a * (1e0 - f);
}

/** @brief Compute linear eccentricity of an ellipse
 *
 * Compute the linear eccentricity of an ellipsoid from the formula
 * \f$ E = \sqrt{a^2 - b^2} \f$
 *
 * @param[in] f Flattening
 * @param[in] a Semi-major axis [me]
 * @return Linear eccentricity \f$E\f$
 */
inline double linear_eccentricity(double a, double f) noexcept {
  const double b = semi_minor(a, f);
  return std::sqrt(a * a - b * b);
}

/** @brief Polar radius of curvature
 *
 * Compute the polar radius of curvature of an ellipsoid (i.e. c) from the
 * formula \f$ c = a^2 / b\f$
 *
 * @param[in] f Flattening
 * @param[in] a Semi-major axis [m]
 * @return Polar radius of curvature
 */
inline constexpr double polar_radius_of_curvature(double a, double f) noexcept {
  const double b = semi_minor(a, f);
  return a * a / b;
}

/** @brief Compute the normal radius of curvature at a given latitude.
 *
 * @param[in] a   The ellipsoid's semi-major axis [m]
 * @param[in] f   Flattening
 * @param[in] lat The latitude [rad]
 * @return Normal radius of curvature at lat [m]
 */
inline double N(double a, double f, double lat) noexcept {
  const double sf = std::sin(lat);
  return a / std::sqrt(1e0 - eccentricity_squared(f) * sf * sf);
}

namespace detail {
/** @brief Compute the normal radius of curvature at a given latitude.
 *
 * @param[in] a   The ellipsoid's semi-major axis [m]
 * @param[in] f   Flattening
 * @param[in] lat The latitude [rad]
 * @param[out] sinlat The computed sin of the given lat (i.e. sin(lat)). Since
 *                we are computing it here, we might as well return it!
 * @return Normal radius of curvature at lat [m]
 */
inline double N(double a, double f, double lat, double &sinlat) noexcept {
  sinlat = std::sin(lat);
  return a / std::sqrt(1e0 - eccentricity_squared(f) * sinlat * sinlat);
}
} // namespace detail

/** @brief Compute the meridional radii of curvature at a given latitude.
 *
 * @param[in] a   The ellipsoid's semi-major axis [m]
 * @param[in] f   Flattening
 * @param[in] lat The latitude in [rad]
 * @return The meridional radius of curvature [m]
 */
inline double M(double a, double f, double lat) noexcept {
  double slat;
  const double Rn = detail::N(a, f, lat, slat);
  return Rn * ((1e0 - eccentricity_squared(f)) /
               (1e0 - eccentricity_squared(f) * slat * slat));
}

/** @brief Compute the geocentric latitude, given a geodetic one for a point
 *        on the ellipsoid (aka, h = 0)
 *
 * The geocentric latitude is the angle between the equatorial plane and the
 * radius from the centre to a point on the surface. The relation between the
 * geocentric latitude(\f$\theta\f$) and the geodetic latitude(\f$\phi\f$) is
 * \f$ \theta (\phi) = tan^{-1} ((1-f)^2 tan(\phi)) \f$
 * The geodetic and geocentric latitudes are equal at the equator and at the
 * poles but at other latitudes they differ by a few minutes of arc.
 * Reference Torge, 2001, Eq. 4.11
 *
 * @param[in] f   The ellipsoid's flattening [-]
 * @param[in] lat The (geodetic) latitude in [rad]
 * @return The geocentric latitude [rad]
 */
inline double geocentric_latitude(double f, double lat) noexcept {
  return std::atan((1e0 - f) * (1e0 - f) * std::tan(lat));
}

/** @brief Compute the parametric or reduced latitude
 *
 * The parametric or reduced latitude, \f$ beta \f$ is defined by the radius
 * drawn from the centre of the ellipsoid to that point Q on the surrounding
 * sphere (of radius a) which is the projection parallel to the Earth's axis
 * of a point P on the ellipsoid at latitude \f$ \phi \f$
 * Reference Torge, 2001, Eq. 4.11
 *
 * @param[in] f   The ellipsoid's flattening
 * @param[in] lat The (geodetic) latitude in radians
 * @return The parametric or reduced latitude at lat in radians
 */
inline double reduced_latitude(double f, double lat) noexcept {
  return std::atan((1e0 - f) * std::tan(lat));
}

} /* namespace core */
} /* namespace dso */

#endif

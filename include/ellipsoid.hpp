/** @file
 * Definition and basic setup of an ellipsoid (kinda) class.
 *
 * This file defines a list of frequently used (referennce) Ellipsoids in
 * geodesy. Each ellipsoid comes with a list of fundamental geometric
 * characteristics (semi-major axis, flattening, name).
 *
 * There are two ways a user can make use of the ellipsoids:
 * * 1. If the ellipsoid of choice is known at compile time, then the
 *      template function/traits can be used. A (per-ellipsoid) specialized
 *      traits class (i.e. dso::ellipsoid_traits) is used to implement the
 *      basic properties of each ellipsoid.
 *      E.g., it the ellipsoid of choice is GRS80, then:
 *      double semi_major = ellipsoid_traits<ellipsoid::grs80>::a;
 *      or
 *      double N = N<ellipsoid::grs80>(lat);
 *      gives the normal radius of curvature at latitude lat.
 *
 * * 2. If the ellipsoid of choice is only kown at runtime, then all relevant
 *      computations/constants can be accessed via the Ellipsoid class; i.e.
 *      Ellipsoid e (ellipsoid::grs80);
 *      double semi_major = e.semi_major();
 *      double N = e.N(lat);
 * Note that the semi-major axis is sometimes reffered to as
 * "equatorial radius" and the semi-minor axis is sometimes reffered to as
 * "polar radius".
 * In the following and when no other qualification is used, latitude is used
 * for Geodetic latitude.
 *
 * [2] Charles F. F. Karney, Algorithms for geodesics, J Geod (2013) 87:43–55
 * [3] https://en.wikipedia.org/wiki/Latitude
 */

#ifndef __DSO_REFERENCE_ELLIPSOID_HPP__
#define __DSO_REFERENCE_ELLIPSOID_HPP__

#include "core/ellipsoid_core.hpp"
#include <type_traits>

namespace dso {

/** @brief A list of well-known reference ellipsoids.
 *
 * For each of these reference ellipsoids, a series of traits (i.e. geometric
 * characteristics) will be specialized later on, using the template class
 * dso::ellipsoid_traits.
 */
enum class ellipsoid : char { grs80, wgs84, pz90 };

/** @brief A class to hold ellipsoid traits (generic case).
 *
 * A (template class) to hold specialized geometric quantities for each
 * of the eumerated elements (i.e. reference ellipsoids) in the
 * dso::ellipsoid enum.
 * I.e., to make any element of dso::ellipsoid usable, specialize this
 * (trait) class.
 *
 * @tparam E  The reference ellipsoid to be specialized (i.e. one of
 *            dso::ellipsoid).
 */
template <ellipsoid E> struct ellipsoid_traits {};

/** @brief A class to hold traits for the GRS-80 (i.e. dso::ellispoid::grs80)
 *  reference ellipsoid.
 *  @see https://en.wikipedia.org/wiki/GRS_80
 */
template <> struct ellipsoid_traits<ellipsoid::grs80> {
  /** Semi-major axis (m). */
  static constexpr double a{6378137e0};
  /** Flattening. */
  static constexpr double f{1e0 / 298.257222101e0};
  /** Reference ellipsoid name. */
  static constexpr const char *n{"GRS80"};
};

/** @brief A class to hold traits for the WGS-84 (i.e. dso::ellispoid::wgs84)
 *  reference ellipsoid.
 *  @see https://en.wikipedia.org/wiki/World_Geodetic_System
 */
template <> struct ellipsoid_traits<ellipsoid::wgs84> {
  /** Semi-major axis (m). */
  static constexpr double a{6378137e0};
  /** Flattening. */
  static constexpr double f{1e0 / 298.257223563e0};
  /** Reference ellipsoid name. */
  static constexpr const char *n{"WGS84"};
};

/** @brief A class to hold traits for the PZ-90 (i.e. dso::ellispoid::pz90)
 * reference ellipsoid.
 * @see
 * http://www.navipedia.net/index.php/Reference_Frames_in_GNSS#GLONASS_reference_frame_PZ-90
 */
template <> struct ellipsoid_traits<ellipsoid::pz90> {
  /** Semi-major axis (m). */
  static constexpr double a{6378136e0};
  /** Flattening. */
  static constexpr double f{1e0 / 298.257839303e0};
  /** Reference ellipsoid name. */
  static constexpr const char *n{"PZ90"};
};

/** @brief Compute the squared eccentricity.
 *
 * @tparam E The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @return Eccentricity squared.
 */
template <ellipsoid E> constexpr double eccentricity_squared() noexcept {
  return core::eccentricity_squared(ellipsoid_traits<E>::f);
}

/** @brief Compute the linear eccentricity.
 *
 * @tparam E  The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @return linear eccentricity
 */
template <ellipsoid E> double linear_eccentricity() noexcept {
  return core::linear_eccentricity(ellipsoid_traits<E>::a,
                                   ellipsoid_traits<E>::f);
}

/** @brief Compute the semi-minor axis (b).
 *
 * @tparam E The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @return Semi-minor axis of the reference ellipsoid in [m].
 */
template <ellipsoid E> constexpr double semi_minor() noexcept {
  return core::semi_minor(ellipsoid_traits<E>::a, ellipsoid_traits<E>::f);
}

/** @brief Compute the polar radius of curvature (c).
 *
 * @tparam E The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @return Polar radius of curvature of the reference ellipsoid in [m].
 */
template <ellipsoid E> constexpr double polar_radius_of_curvature() noexcept {
  return core::polar_radius_of_curvature(ellipsoid_traits<E>::a,
                                         ellipsoid_traits<E>::f);
}

/** @brief Compute the third flattening
 *
 * @tparam E The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @return Third flattening
 */
template <ellipsoid E> constexpr double third_flattening() noexcept {
  return core::third_flattening(ellipsoid_traits<E>::f);
}

/** @brief Compute the normal radius of curvature at a given latitude.
 *
 * @tparam E  The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @param[in] lat The latitude in [rad].
 * @return    The normal radius of curvature in [m].
 */
template <ellipsoid E> double N(double lat) noexcept {
  return core::N(ellipsoid_traits<E>::a, ellipsoid_traits<E>::f, lat);
}

/** @brief Compute the geocentric latitude at some geodetic latitude on the
 *        ellipsoid (that is ellipsoidal height = 0).
 *
 * @tparam  E The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @param[in] lat The geodetic latitude [rad]
 * @return    The geocentric latitude [rad]
 *
 * @note If the point has height != 0 (aka, is not ON the ellipsoid), the use
 *  the version whilch also takes height as input.
 */
template <ellipsoid E> double geocentric_latitude(double lat) noexcept {
  return core::geocentric_latitude(ellipsoid_traits<E>::f, lat);
}

/** @brief Compute the geocentric latitude at some geodetic latitude, for a
 *        point with a given ellipsoidal height.
 *
 * @tparam E  The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @param[in] lat The geodetic latitude [rad]
 * @param[in] h   The ellipsoidal height [m]
 * @return The geocentric latitude [rad]
 *
 * @note If the point has height = 0 (aka, is ON the ellipsoid), then use
 *  the version whilch does not take height as input.
 *
 * @see
 * https://www.mathworks.com/help/aeroblks/geodetictogeocentriclatitude.html
 *      note that in the polar axis distance there is a sin when it should
 *      have been a cos.
 */
template <ellipsoid E>
double geocentric_latitude(double lat, double h) noexcept {
  constexpr const double ecc2 = eccentricity_squared<E>();
  const double slat = std::sin(lat);
  const double clat = std::cos(lat);
  const double Rn = N<E>(lat);
  const double rho = (Rn + h) * clat;
  const double z = (Rn * (1e0 - ecc2) + h) * slat;
  return std::atan(z / rho);
}

/** @brief Compute the parametric or reduced latitude at some geodetic latitude
 *        on the ellipsoid
 *
 * @tparam  E The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @param[in] lat The geodetic latitude in [rad]
 * @return    The parametric or reduced latitude in [rad]
 */
template <ellipsoid E> double reduced_latitude(double lat) noexcept {
  return core::reduced_latitude(ellipsoid_traits<E>::f, lat);
}

/** @brief Compute the meridional radii of curvature at a given latitude
 *
 * @tparam  E The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @param[in] lat The latitude in [rad].
 * @return The meridional radius of curvature in [m].
 */
template <ellipsoid E> double M(double lat) noexcept {
  return core::M(ellipsoid_traits<E>::a, ellipsoid_traits<E>::f, lat);
}

/** @brief Compute mean earth radius; in geophysics, the International Union of
 *        Geodesy and Geophysics (IUGG) defines the mean radius (denoted R1)
 *        to be: \f$ R1 = (2\alpha + b) / 3 \f$
 *
 * @tparam E The reference ellipsoid (i.e. one of dso::ellipsoid).
 * @return The mean earth radius (as defined by IUGG for ellipsoid E) in [m]
 *
 * @see https://en.wikipedia.org/wiki/Earth_radius#Mean_radius
 */
template <ellipsoid E> constexpr double mean_earth_radius() noexcept {
  return 2e0 * ellipsoid_traits<E>::a / 3e0 + semi_minor<E>() / 3e0;
}

/** @brief Arc length of an infinitesimal element on the meridian
 *
 * @warning This formula is valid for infinitesimal latitude differences
 * @see https://en.wikipedia.org/wiki/Meridian_arc
 *
 * @param[in] lat Latitude of the reference point in [rad]
 * @param[in] dlat Lattitude difference, i.e. arc length on the meridian,
 *            in [rad]
 * @return  Arc length (on meridian) in [m]
 */
template <ellipsoid E>
double infinitesimal_meridian_arc(double lat, double dlat) noexcept {
  const double Rm = M<E>(lat);
  return Rm * dlat;
}

/** @brief Arc length on parallel
 *
 * @param[in] lat  Latitude of the parallel [rad]
 * @param[in] dlon Longtitude difference [rad]
 * @return Arc length (on parallel) in [m]
 */
template <ellipsoid E>
double parallel_arc_length(double lat, double dlon) noexcept {
  const double Rn = N<E>(lat);
  const double cosf = std::cos(lat);
  return Rn * cosf * dlon;
}

/** @class Ellipsoid
 *
 * A class to represent a reference ellipsoid. An ellipsoid is defined by
 * two parameters, namely:
 * * semi-major axis, \f$ \alpha \f$ aka the equatorial radius of the ellipsoid
 * * flattening, f aka \f$ f = \frac{\alpha - \beta}{\alpha} \f$
 *
 * Users can construct the commonly used ellipsoids in Geodesy (grs80,
 * wgs84 and pz90) via the dso::ellipsoid enums, or any other ellipsoid
 * of choice, by passing in the fundamental arguments (a and f).
 */
class Ellipsoid {
public:
  /** @brief  Constructor from an dso::ellipsoid enum
   * @param[in] e An dso::ellipsoid; fundamental geometric constants are
   *   automatically assigned via the dso::ellipsoid_traits class.
   */
  explicit constexpr Ellipsoid(ellipsoid e) noexcept : __a(0e0), __f(0e0) {
    switch (e) {
    case ellipsoid::grs80:
      __a = ellipsoid_traits<ellipsoid::grs80>::a;
      __f = ellipsoid_traits<ellipsoid::grs80>::f;
      break;
    case ellipsoid::wgs84:
      __a = ellipsoid_traits<ellipsoid::wgs84>::a;
      __f = ellipsoid_traits<ellipsoid::wgs84>::f;
      break;
    case ellipsoid::pz90:
      __a = ellipsoid_traits<ellipsoid::pz90>::a;
      __f = ellipsoid_traits<ellipsoid::pz90>::f;
      break;
    }
  }

  /** @brief  User-defined instance.
   * @param[in] a The semi-major axis [m]
   * @param[in] f The flattening
   */
  constexpr Ellipsoid(double a, double f) noexcept : __a(a), __f(f) {};

  /** @brief et the semi-major axis \f$ \alpha \f$ */
  constexpr double semi_major() const noexcept { return __a; }

  /** @brief  Get the flatteninga */
  constexpr double flattening() const noexcept { return __f; }

  /** @brief  Get the squared eccentricity \f$ e^2 \f$ */
  constexpr double eccentricity_squared() const noexcept {
    return core::eccentricity_squared(__f);
  }

  /** @brief  Get the semi-minor axis \f$ \beta \f$ */
  constexpr double semi_minor() const noexcept {
    return core::semi_minor(__a, __f);
  }

  /** @brief Get the third flattening \f$ n \f$ */
  constexpr double third_flattening() const noexcept {
    return core::third_flattening(__f);
  }

  /** @brief Compute the geocentric latitude at some (geodetic) latitude */
  double geocentric_latitude(double lat) const noexcept {
    return core::geocentric_latitude(__f, lat);
  }

  /** @brief Compute the reduced latitude at some (geodetic) latitude */
  double reduced_latitude(double lat) const noexcept {
    return core::reduced_latitude(__f, lat);
  }

  /** @brief  Compute the normal radius of curvature at a given latitude */
  double N(double lat) const noexcept { return core::N(__a, __f, lat); }

  /** @brief  Compute the meridional radii of curvature at a given latitude */
  double M(double lat) const noexcept { return core::M(__a, __f, lat); }

private:
  double __a, __f;
}; /* class Ellipsoid */

} /* namespace dso */

#endif

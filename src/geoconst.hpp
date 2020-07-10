///
/// @file geoconst.hpp
///
/// @brief Fundamental constants frequently used within geodetic calculations.
///

#ifndef __NGPT_GEOCONST_HPP__
#define __NGPT_GEOCONST_HPP__

#include <cmath>

namespace ngpt {
/// The value of pi.
#if defined(__GNUC__) && !defined(__llvm__)
constexpr double DPI{std::atan(1e0) * 4e0};
#else
#define _USE_MATH_DEFINES
constexpr double DPI{M_PI};
#endif

/// The value of 2 * pi.
constexpr double D2PI{2e0 * DPI};

/// Degrees to Radians coefficient.
constexpr double DEG2RAD{DPI / 180e0};

/// Radians to Degrees coefficient.
constexpr double RAD2DEG{180e0 / DPI};

/// mas to radians factor aka θrad = θmas*MAS2RAD
constexpr double MAS2RAD{4.847309743e-9};

} // namespace ngpt

#endif

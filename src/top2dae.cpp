#include "geodesy.hpp"
#include "units.hpp"

double dso::top2dae(const Eigen::Matrix<double, 3, 1> &enu, double &azimouth,
                    double &elevation) {
  const double e = enu(0);
  const double n = enu(1);
  const double u = enu(2);

  const double rho2 = e * e + n * n;
  const double rho = std::sqrt(rho2);

  // azimouth in [0,2π]
  azimouth = norm_angle<dso::AngleUnit::Radians>(std::atan2(e, n));

  // elevation angle [0-π]
  elevation = std::atan(u / rho);

  // distance
  return enu.norm();
}

double dso::top2dae(const Eigen::Matrix<double, 3, 1> &enu, double &azimouth,
                    double &elevation, Eigen::Matrix<double, 3, 1> &dAdr,
                    Eigen::Matrix<double, 3, 1> &dEdr) {
  const double e = enu(0);
  const double n = enu(1);
  const double u = enu(2);

  const double rho2 = e * e + n * n;
  const double rho = std::sqrt(rho2);

  // azimouth in [0,2π]
  azimouth = norm_angle<dso::AngleUnit::Radians>(std::atan2(e, n));

  // elevation angle [0-π]
  elevation = std::atan(u / rho);

  // partials
  dAdr(0) = n / rho2;
  dAdr(1) = -e / rho2;
  dAdr(2) = 0e0;

  const double R2 = rho2 + u * u;
  dEdr(0) = -e * u / rho / R2;
  dEdr(1) = -n * u / rho / R2;
  dEdr(2) = rho / R2;

  // distance
  return enu.norm();
}

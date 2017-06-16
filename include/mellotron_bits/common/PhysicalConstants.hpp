#ifndef PHYSICAL_CONSTANTS_HPP
#define PHYSICAL_CONSTANTS_HPP

#include <boost/math/constants/constants.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/atomic-nuclear_constants.hpp>

namespace mellotron {

/*!
 *  \file   PhysicalConstants.hpp
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2016-12-08
 *  \brief  Easy access to useful mathematical/physical constants.
 */

namespace constants {

namespace math {

  const double pi = boost::math::constants::pi<double>();                                                                                ///< The circle constant.

} // namespace math

namespace physics {

  const double electron_mass      = boost::units::si::constants::codata::m_e / boost::units::si::kilogram;                               ///< Electron rest mass in kg.
  const double electron_charge    = boost::units::si::constants::codata::e / boost::units::si::coulomb;                                  ///< Elementary charge in C.
  const double c                  = boost::units::si::constants::codata::c / boost::units::si::meter * boost::units::si::second;         ///< Speed of light in m/s.
  const double epsilon_0          = boost::units::si::constants::codata::epsilon_0 / boost::units::si::farad * boost::units::si::meter;  ///< Vacuum permittivity in F/m.
  const double hbar               = boost::units::si::constants::codata::hbar / boost::units::si::joule / boost::units::si::second;      ///< hbar.
  const double alpha              = boost::units::si::constants::codata::alpha / boost::units::si::dimensionless();                      ///< Fine structure constant.
  const double UNIT_TIME_QED      = 1.28808867e-21;                                                                                      ///< Unit time in QED units, hbar/(m_e*c^2).
  const double UNIT_TIME_EV       = 6.582119e-16;                                                                                        ///< Unit time in eV units.

} // namespace physics

} // namespace constants

} // namespace mellotron

#endif // PHYSICAL_CONSTANTS_HPP
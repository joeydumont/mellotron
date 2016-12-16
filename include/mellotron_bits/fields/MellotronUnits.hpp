#ifndef MELLOTRON_UNITS_HPP
#define MELLOTRON_UNITS_HPP

#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>

#include "mellotron_bits/common/PhysicalConstants.hpp"


namespace mellotron {

/*!
 *  \class  MellotronUnits
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2016-12-07
 *  \brief  Defines the base quantities in our unit system.
 *
 * In the Mellotron, we use a system of units we dub electronic units.
 * In short, the fields are normalized via E' = |e|/(m_e omega_0 c) E
 * where e is the elementary charge, m_e the mass of the electron, c the speed
 * of light and omega_0 a characteristic frequency of the field. Since omega_0
 * can change, in contradistinction to the electron mass, our unit system varies
 * with the frequency. This class computes the unit_length, unit_momentum, unit_energy
 * and other such quantities depending on the value of omega_0.
 */

typedef enum {si,qed,ev} unit_system;

class MellotronUnits
{
public:
  /// Defines the unit system.
  MellotronUnits(double my_omega_0, unit_system units = si)
  {
    // We compute omega_0 in SI units.
    switch (units)
    {
      case (si) : {omega_0_SI = my_omega_0;break;}
      case (qed): {omega_0_SI = my_omega_0/constants::physics::UNIT_TIME_QED;break;}
      case (ev) : {omega_0_SI = my_omega_0/constants::physics::UNIT_TIME_EV;break;}
      default:    {std::cout << "Unknown units." << std::endl;throw;break;}
    }

    UNIT_LENGTH   = constants::physics::c/omega_0_SI;
    UNIT_MOMENTUM = constants::physics::electron_mass*constants::physics::c;
    UNIT_TIME     = 1.0/omega_0_SI;
    UNIT_ENERGY   = constants::physics::epsilon_0*std::pow(constants::physics::electron_mass,2)*std::pow(constants::physics::c,5)/(std::pow(constants::physics::electron_charge,2)*omega_0_SI);
  }

  // Variable data.
  double omega_0_SI;       ///< Characteristic frequency of the beam, in SI units.
  double UNIT_LENGTH;      ///< Unit length in Mellotron units.
  double UNIT_MOMENTUM;    ///< Unit momentum in Mellotron units.
  double UNIT_TIME;        ///< Unit time in Mellotron units.
  double UNIT_ENERGY;      ///< Unit energy in Mellotron units.
};

} // namespace mellotron

#endif // MELLOTRON_UNITS_HPP

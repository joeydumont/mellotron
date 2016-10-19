#ifndef PARTICLE_BONES_HPP
#define PARTICLE_BONES_HPP

#include <armadillo>
#include <cmath>

/*!
 *  \class  Particle
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2016-09-30
 *  \brief  Defines the Lorentz force equation acting on a particle with given
 *         properties.
 *
 * This class defines the properties of a charged particle. We also define
 * its equation of motion, i.e. the Lorentz equation plus, possibly, radiation
 * reaction terms.
 */
template <class FieldModel>
class Particle
{
public:

  /// Constructor sets the physical properties of the particle.
  Particle(const double my_charge, const double my_mass, FieldModel& my_field_model, const std::string my_radiation_reaction = std::string("no_rr"));

  /// Computation of the field tensor at a given point in spacetime.
  void ComputeFieldTensor(const double t, const double x, const double y, const double z);

  // Accessor functions of the electromagnetic fields.
  arma::colvec::fixed<3> GetElectricField(){return electric_field;} ///< Returns the stored electric field.
  arma::colvec::fixed<3> GetMagneticField(){return magnetic_field;} ///< Returns the stored magnetic field.

  // Accessor functions of the particle parameters.
  double GetChi(){return chi;}                                      ///< Returns the dynamical quantum parameter of the particle.

  /// Overloading of the () operator for use with Boost.odeint.
  void operator()(const arma::colvec::fixed<8>& x, arma::colvec::fixed<8> &dxdt, const double t);

protected:

        FieldModel               &  field_model;           ///< Object that contains a ComputeFieldComponents routine.

        arma::colvec::fixed<3>      electric_field;        ///< Electric field at a given point in spacetime.
        arma::colvec::fixed<3>      magnetic_field;        ///< Magnetic field at a given point in spacetime.

  const double                      charge;                ///< Charge of the particle, in QED units.
  const double                      mass;                  ///< Mass of the particle, in QED units.

  const std::string                 radiation_reaction;    ///< Determines the model of radiation reaction we employ, if at all.

  const double                      alpha;                 ///< Fine-structure constant.

        double                      chi_sq;                ///< Lorentz invariant along the trajectory, squared.
        double                      chi;                   ///< Lorentz invariant along the trajectory.
};

#endif // PARTICLE_BONES_HPP

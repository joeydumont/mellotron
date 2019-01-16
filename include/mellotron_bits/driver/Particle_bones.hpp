#ifndef PARTICLE_BONES_HPP
#define PARTICLE_BONES_HPP

#include <armadillo>
#include <cmath>

#include "mellotron_bits/driver/Envelope_bones.hpp"

namespace mellotron {

/// enum to choose the radiation reaction model.
enum RadiationReactionModel {NoRR, LandauLifshitz, LandauLifshitzQuantumCorrection};

/*!
 *  \class  Particle
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2016-09-30
 *  \brief  Defines the Lorentz force equation acting on a particle with given
 *          properties.
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
  /// RR makes it so the unit system must be passed to the Particle.
  Particle(const double                    my_charge,
           const double                    my_mass,
                 FieldModel             &  my_field_model,
                 MellotronUnits         &  my_units,
                 Envelope               &  my_envelope,
           const RadiationReactionModel    my_radiation_reaction = NoRR);

  /// Computation of the field tensor at a given point in space-time.
  void ComputeFieldTensor(const double t, const double x, const double y, const double z);

  // Accessor functions of the electromagnetic fields.
  arma::colvec::fixed<3> GetElectricField(){return electric_field;} ///< Returns the stored electric field.
  arma::colvec::fixed<3> GetMagneticField(){return magnetic_field;} ///< Returns the stored magnetic field.

  // Accessor functions of the particle parameters.
  double GetCharge(){return charge;}                                ///< Returns the charge of the particle.
  double GetMass(){return mass;}                                    ///< Returns the mass of the particle.
  double GetChi(){return chi;}                                      ///< Returns the dynamical quantum parameter of the particle.

  /// Accessor function of the unit system.
  MellotronUnits & GetUnitSystem(){return unit_system;}

  // Utility function to set the initial conditions.
  void SetInitConditions(arma::colvec::fixed<8>& x, double x_init, double y_init, double z_init, double px_init, double py_init, double pz_init, double t_init);

  /// Overloading of the () operator for use with Boost.odeint.
  void operator()(const arma::colvec::fixed<8>& x, arma::colvec::fixed<8> &dxdt, const double t);

  /// Function that computes the Lorentz force.
  void ComputeLorentzForce(const arma::colvec::fixed<8> & x, arma::colvec::fixed<8> & dxdt, const double t);

protected:

  const double                      charge;                ///< Charge of the particle, multiple of the elementary charge.
  const double                      mass;                  ///< Mass of the particle, multiple of the electron mass.

        FieldModel               &  field_model;           ///< Object that contains a ComputeFieldComponents routine.

        arma::colvec::fixed<3>      electric_field;        ///< Electric field at a given point in spacetime.
        arma::colvec::fixed<3>      magnetic_field;        ///< Magnetic field at a given point in spacetime.


        MellotronUnits           &  unit_system;           ///< Contains information about the unit system used. Useful for RR.
        Envelope                 &  envelope;              ///< Object that contains the electromagnetic field envelope
  const RadiationReactionModel      radiation_reaction;    ///< Determines the model of radiation reaction we employ, if at all.

        double                      chi_sq;                ///< Lorentz invariant along the trajectory, squared.
        double                      chi;                   ///< Lorentz invariant along the trajectory.

};

/*!
 *  \class  ParticleIonized
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2017-06-20
 *  \brief  Defines the Lorentz force equation acting on a particle with given
 *          properties. Includes a field threshold below which to Lorentz force
 *          is applied.
 *
 * This class defines the properties of a charged particle. We also define
 * its equation of motion, i.e. the Lorentz equation plus, possibly, radiation
 * reaction terms. A field threshold below which the Lorentz force is assumed
 * is vanish is also implemented. This simulates, in some approximate way,
 * the ionization process.
 */
template <class FieldModel>
class ParticleIonized : public Particle<FieldModel>
{
public:

  /// Constructor sets the physical properties of the particle.
  /// It also sets the field threshold above which we apply the Lorentz force.
  /// To have no threshold, use Particle, of a negative value of field_threshold.
  /// RR makes it so the unit system must be passed to the Particle.
  ParticleIonized(const double                     my_charge,
                  const double                     my_mass,
                        FieldModel              &  my_field_model,
                        MellotronUnits          &  my_units,
                        double                     my_field_threshold,
                        Envelope                &  my_envelope,
                  const RadiationReactionModel     my_radiation_reaction = NoRR);

  /// Accessor function to check whether the particle has been "ionized".
  bool GetIonized(){return ApplyLorentzForce;}

  /// Overloading of the () operator for use with Boost.odeint.
  void operator()(const arma::colvec::fixed<8>& x, arma::colvec::fixed<8> &dxdt, const double t);

protected:
        double                      field_threshold;       ///< Field strength above which we apply the Lorentz force.
        bool                        ApplyLorentzForce;     ///< Flag to determine if the particle has been ionized yet, and to apply the Lorentz force if so.

};

} // namespace mellotron

#endif // PARTICLE_BONES_HPP

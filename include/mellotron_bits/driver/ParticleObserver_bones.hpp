#ifndef PARTICLE_OBSERVER_BONES_HPP
#define PARTICLE_OBSERVER_BONES_HPP

#include <armadillo>
#include <boost/functional/hash.hpp>
#include <cmath>
#include <hdf5.h>

namespace mellotron {

/*!
 *  \class  ParticleObserver
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2016-09-30
 *  \brief  Observes and outputs the trajectory of a particle in an electromagnetic field.
 *
 * Defines the output of Mellotron. This class observes a single particle as it moves through
 * space. Contains a reference to a Particle object as to query the value of the electric and
 * magnetic field.
 */
template <class FieldModel>
class ParticleObserver
{
public:

  /// Sets the particle to follow.
  ParticleObserver(Particle<FieldModel>& my_particle);

  /// Overloading of the () operator for use with Boost.odeint.
  void operator()(const arma::colvec::fixed<8>& x, double t);

  /// Outputs in a HDF5 file.
  void OutputData();



       Particle<FieldModel>  &  particle;              ///< Particle object that moves through spacetime.

       arma::mat                position;              ///< Record of the particle's position.
       arma::mat                momentum;              ///< Record of the particle's momentum.
       std::vector<double>      gamma;                 ///< Record of the particle's gamma factor.
       std::vector<double>      chi;                   ///< Record of the particle's Lorentz invariant.
       std::vector<double>      times;                 ///< Record of the time values.

       arma::mat                electric_field;        ///< Record of the electric field along the particle's trajectory.
       arma::mat                magnetic_field;        ///< Record of the magnetic field along the particle's trajectory.
};

/*!
 *  \class  ParticleObserverIonized
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2017-06-21
 *  \brief  Observes and outputs the trajectory of an "ionized" particle in an electromagnetic field.
 *
 * This class observes a single particle as it moves through
 * space. Contains a reference to a ParticleIonized object as to query the value of the electric and
 * magnetic field. It only outputs the data if the particle has been "ionized" during the simulation,
 * i.e. that the field surpassed the threshold at some point during the simulation.
 */
template <class FieldModel>
class ParticleObserverIonized : public ParticleObserver<FieldModel>
{
public:

  /// Sets the particle to follow.
  ParticleObserverIonized(ParticleIonized<FieldModel>& my_particle);

  /// Outputs in a HDF5 file.
  void OutputData();

       ParticleIonized<FieldModel>  &  particle_ion;              ///< Particle object that moves through spacetime.
};

} // namespace mellotron

#endif // PARTICLE_OBSERVER_BONES_H

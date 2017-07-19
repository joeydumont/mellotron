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
  /// Also sets an initial size for the data structures that hold information
  /// about the particle.
  ParticleObserver(Particle<FieldModel>& my_particle, const int my_init_size = 100);

  /// Overloading of the () operator for use with Boost.odeint.
  void operator()(const arma::colvec::fixed<8>& x, double t);

  /// Outputs in a pair of HDF5 and XDMF files.
  void OutputData();

  /// Creates the HDF5 and outputs the data.
  void  GenerateHDF5();

        Particle<FieldModel>  &  particle;              ///< Particle object that moves through spacetime.

  const int                      init_size;             ///< Initial size of the data structures. Removes some unnecessary arma::resize() calls.

        arma::mat                position;              ///< Record of the particle's position.
        arma::mat                momentum;              ///< Record of the particle's momentum.
        std::vector<double>      gamma;                 ///< Record of the particle's gamma factor.
        std::vector<double>      chi;                   ///< Record of the particle's Lorentz invariant.
        std::vector<double>      times;                 ///< Record of the time values.

        arma::mat                electric_field;        ///< Record of the electric field along the particle's trajectory.
        arma::mat                magnetic_field;        ///< Record of the magnetic field along the particle's trajectory.

        uint                     step_counter;          ///< Counts the number of steps that were taken.
};

} // namespace mellotron

#endif // PARTICLE_OBSERVER_BONES_H

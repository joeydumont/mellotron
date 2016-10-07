#ifndef PARTICLE_OBSERVER_H
#define PARTICLE_OBSERVER_H

#include <armadillo>
#include <cmath>

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

  /// Outputs in a pair of HDF5 and XDMF files.
  void OutputData();

protected:

        Particle<FieldModel>  &  particle;              ///< Particle object that moves through spacetime.

        arma::mat                position;              ///< Record of the particle's position.
        arma::mat                momentum;              ///< Record of the particle's momentum.
        std::vector<double>      gamma;                 ///< Record of the particle's gamma factor.
        std::vector<double>      chi;                   ///< Record of the particle's Lorentz invariant.

        arma::mat                electric_field;        ///< Record of the electric field along the particle's trajectory.
        arma::mat                magnetic_field;        ///< Record of the magnetic field along the particle's trajectory.
};

#endif // PARTICLE_OBSERVER_H

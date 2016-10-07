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
  ParticleObserver(Particle& my_particle)
  : particle(my_particle)
  {}

  /// Overloading of the () operator for use with Boost.odeint.
  void operator()(const arma::colvec::fixed<8>& x, double t);

protected:

        Particle                &  particle;              ///< Particle object that moves through spacetime.

        arma::colvec::fixed<3>     electric_field;        ///< Electric field at a given point in spacetime.
        arma::colvec::fixed<3>     magnetic_field;        ///< Magnetic field at a given point in spacetime.
};

#endif // PARTICLE_H

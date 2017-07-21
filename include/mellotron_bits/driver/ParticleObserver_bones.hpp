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
  ParticleObserver(Particle<FieldModel>& my_particle, const uint my_init_size = 100);

  /// Overloading of the () operator for use with Boost.odeint.
  void operator()(const arma::colvec::fixed<8>& x, double t);

  /// Outputs in a pair of HDF5 and XDMF files.
  void OutputData();

  /// Writes all the data in a given HDF5 group. Useful to append other write functions in derived classes.
  void WriteAllData(hid_t group_id);

  /// Writes the temporal data points.
  void WriteTimes(hid_t group_id);

  /// Writes gamma.
  void WriteGamma(hid_t group_id);

  /// Writes chi.
  void WriteChi(hid_t group_id);

  /// Writes the state vector (position and momentum).
  void WriteStateVector(hid_t group_id);

  /// Write the electric and magnetic field.
  void WriteElectromagneticField(hid_t group_id);

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

protected:

  /// Utility function that sets the right properties for HDF5 output.
  void SetHDF5Properties(hid_t & dataspace_id, hid_t & plist_id, const int dim, const hsize_t * size, const hsize_t * chunk_size, const uint compression_level);

};

/*!
 *  \class  ParticleObserverLienardWiechert
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2016-09-30
 *  \brief  Observes and outputs the trajectory of a particle in an electromagnetic field.
 *
 * Adds the fields created by the moving charges computed via the Liénard-Wiechert potentials
 * to ParticleObserver.
 */
template <class FieldModel>
class ParticleObserverLienardWiechert : public ParticleObserver<FieldModel>
{
public:

  /// Sets the particle to follow.
  /// We must also declare the radius of the sphere on which
  /// the Liénard-Wiechert field are computed, and the number of points
  /// in theta and phi.
  ParticleObserverLienardWiechert(       Particle<FieldModel>  &  my_particle,
                                  const  double                   my_radius              = 1e4,
                                  const  unsigned int             my_number_points_theta = 50,
                                  const  unsigned int             my_number_points_phi   = 50,
                                  const  unsigned int             my_init_size           = 100);

  /// We redefine the () operator to append the calculation of the emitted radiation.
  void operator()(const arma::colvec::fixed<8>& x, double t);


  const double                   radius;                ///< Radius of the detection sphere.
  const unsigned int             number_points_theta;   ///< Number of discrete points in theta.
  const unsigned int             number_points_phi;     ///< Number of discrete points in phi.
};

} // namespace mellotron

#endif // PARTICLE_OBSERVER_BONES_H

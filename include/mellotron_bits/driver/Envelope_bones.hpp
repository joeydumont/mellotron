#ifndef ENVELOPE_BONES_HPP
#define ENVELOPE_BONES_HPP

#include <cmath>

#include "mellotron_bits/common/PhysicalConstants.hpp"

namespace mellotron {

/*!
 *  \class  Envelope
 *  \author Francois Fillion-Gourdeau      <francois.fillion@emt.inrs.ca>
 *  \since  2018-12-21
 *  \brief  Defines en envelope over the time dependence of the electromagnetic
 *          field
 *
 * This class defines a time-dependent envelope with value in [0,1], that is a
 * multiplicative function that multiplies the  electromagnetic field. It it used
 * mostly to guarantee that the field is zero at the initial time. Hopefully,
 * the implemented envelope models have negligible effect on the spectral
 * components of the electromagnetic field.
 */

class Envelope
{
public:

  /// Constructor
  Envelope();
  /// Destructor
  ~Envelope();

  /// Function that returns the value of the envelope at time t
  virtual double Value(double my_t) = 0;


protected:


};

/*!
 *  \class  NoEnvelope
 *  \author Francois Fillion-Gourdeau      <francois.fillion@emt.inrs.ca>
 *  \since  2018-12-21
 *  \brief  No envelope, i.e. the value is always one.
 */

class NoEnvelope : public Envelope
{
public:

  /// Constructor
  NoEnvelope();
  /// Destructor
  ~NoEnvelope();

  /// Function that returns the value of the envelope at time t
  double Value(double my_t);

protected:

};

/*!
 *  \class  EnvelopeHann
 *  \author Francois Fillion-Gourdeau      <francois.fillion@emt.inrs.ca>
 *  \since  2018-12-21
 *  \brief  Defines an envelope which is a Hann function: g(t) = cos^2(Omega(t-t0))
 */

class EnvelopeHann : public Envelope
{
public:

  /// Constructor that sets the parameter
  EnvelopeHann(double my_initial_time,
               double my_total_time);
  /// Destructor
  ~EnvelopeHann();

  /// Function that returns the value of the envelope at time t
  double Value(double my_t);

protected:

  double initial_time;    ///< Initial time of the simulation
  double total_time;      ///< Total time of the simulation
  double omega;           ///< Frequency of the cos function

};



} // namespace mellotron

#endif // ENVELOPE_BONES_HPP

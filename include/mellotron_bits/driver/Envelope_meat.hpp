#ifndef ENVELOPE_MEAT_HPP
#define ENVELOPE_MEAT_HPP

namespace mellotron {


inline
Envelope::Envelope(){}

inline
Envelope::~Envelope(){}


inline
NoEnvelope::NoEnvelope(){}

inline
NoEnvelope::~NoEnvelope(){}

inline
double NoEnvelope::Value(double my_t){return 1.0;}


inline
EnvelopeHann::EnvelopeHann(double my_initial_time,
                           double my_total_time)
: initial_time(my_initial_time)
, total_time(my_total_time)
{
  omega = constants::math::pi/total_time;
}

inline
EnvelopeHann::~EnvelopeHann(){}

inline
double EnvelopeHann::Value(double my_t)
{
  double cos_value = std::cos(omega*my_t);
  return cos_value*cos_value;
}



} // namespace mellotron

#endif // ENVELOPE_MEAT_HPP

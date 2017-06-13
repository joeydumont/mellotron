#ifndef PARTICLE_MEAT_HPP
#define PARTICLE_MEAT_HPP

namespace mellotron {

template <class FieldModel>
inline
Particle<FieldModel>::Particle(const double          my_charge,
                               const double          my_mass,
                               FieldModel         &  my_field_model,
                               const std::string     my_radiation_reaction)
: charge(my_charge)
, mass(my_mass)
, field_model(my_field_model)
, radiation_reaction(my_radiation_reaction)
, alpha(7.2973525664e-3)
{}

template <class FieldModel>
inline
void
Particle<FieldModel>::ComputeFieldTensor(const double t,
                                         const double x,
                                         const double y,
                                         const double z)
{
  std::array<double,6> field = this->field_model.ComputeFieldComponents(t,x,y,z);

  for (uint i=0; i<3; i++)
  {
    this->electric_field(i) = field[i];
    this->magnetic_field(i) = field[i+3];
  }
}

template <class FieldModel>
inline
void
Particle<FieldModel>::SetInitConditions(arma::colvec::fixed<8>  &x,
                                        const double            x_init,
                                        const double            y_init,
                                        const double            z_init,
                                        const double            px_init,
                                        const double            py_init,
                                        const double            pz_init,
                                        const double            t_init)
{
  // We compute the initial gamma factor.
  double p_init_sq = std::pow(px_init,2)+std::pow(py_init,2)+std::pow(pz_init,2);
  double gamma     = std::sqrt(1.0+p_init_sq/std::pow(mass,2));

  // We set the values in the vector.
  x[0] = t_init;
  x[1] = x_init;
  x[2] = y_init;
  x[3] = z_init;
  x[4] = mass*gamma;
  x[5] = px_init;
  x[6] = py_init;
  x[7] = pz_init;
}

template <class FieldModel>
inline
void
Particle<FieldModel>::operator() (const arma::colvec::fixed<8> &x,
                                        arma::colvec::fixed<8> &dxdt,
                                  const double t)
{
  // We define some auxiliary variables.
  arma::colvec position = x.subvec(1,3);
  arma::colvec momentum = x.subvec(5,7);

  // We compute the field tensor.
  ComputeFieldTensor(x[0],x[1],x[2],x[3]);

  // We compute the Lorentz gamma factor.
  double gamma              = std::sqrt(1.0+std::pow(arma::norm(momentum,2)/mass,2));
  double rel_mass           = gamma*mass;
  double charge_to_rel_mass = charge/rel_mass;

  // Set the position differentials.
  dxdt.subvec(0,3) = x.subvec(4,7)/rel_mass;

  // Compute the force.
  double pdotE         = arma::dot(momentum,electric_field);
  arma::colvec lorentz = gamma*electric_field + arma::cross(momentum,magnetic_field);
  chi_sq               = std::pow(arma::norm(lorentz,2),2)-std::pow(pdotE,2);
  chi                  = std::sqrt(chi_sq);

  // Set the momentum differentials.
  dxdt(4)              = charge_to_rel_mass*arma::dot(momentum,electric_field);
  dxdt.subvec(5,7)     = charge*electric_field+charge_to_rel_mass*arma::cross(momentum,magnetic_field);

  // Radiation reaction effects. TODO: CORRECT CHI VALUE.
  if (radiation_reaction == std::string("ll_first_term"))
  {
    dxdt.subvec(4,7) -= 2.0*alpha*std::pow(charge,4)/(3.0*gamma*std::pow(mass,5))*chi_sq*x.subvec(4,7);
  }

  if (radiation_reaction == std::string("ll_first_term_quantum_correction"))
  {
    double quantum_factor = std::pow(1+18.0*chi+69.0*chi_sq*73.0*std::pow(chi,3)+5.806*std::pow(chi_sq,2),-1.0/3.0);
    dxdt.subvec(4,7) -= quantum_factor*2.0*alpha*std::pow(charge,4)/(3.0*gamma*std::pow(mass,5))*chi_sq*x.subvec(4,7);
  }
}

} // namespace mellotron

#endif // PARTICLE_MEAT_HPP
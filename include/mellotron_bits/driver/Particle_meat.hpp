#ifndef PARTICLE_MEAT_HPP
#define PARTICLE_MEAT_HPP

namespace mellotron {

template <class FieldModel>
inline
Particle<FieldModel>::Particle(const double                     my_charge,
                               const double                     my_mass,
                                     FieldModel              &  my_field_model,
                                     MellotronUnits          &  my_units,
                               const RadiationReactionModel     my_radiation_reaction)
: charge(my_charge)
, mass(my_mass)
, field_model(my_field_model)
, unit_system(my_units)
, radiation_reaction(my_radiation_reaction)
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
  // We compute the field tensor.
  ComputeFieldTensor(x[0],x[1],x[2],x[3]);

  ComputeLorentzForce(x, dxdt, t);
}

template <class FieldModel>
inline
void
Particle<FieldModel>::ComputeLorentzForce(const arma::colvec::fixed<8> &x,
                                                arma::colvec::fixed<8> &dxdt,
                                          const double t)
{
  // We define some auxiliary variables.
  arma::colvec position     = x.subvec(1,3);
  arma::colvec momentum     = x.subvec(5,7);

  // We compute the Lorentz gamma factor.
  double gamma              = std::sqrt(1.0+std::pow(arma::norm(momentum,2)/mass,2));
  double rel_mass           = gamma*mass;
  double charge_to_rel_mass = charge/rel_mass;

  // Set the position differentials.
  dxdt.subvec(0,3)          = x.subvec(4,7)/rel_mass;

  // Compute the force.
  double pdotE              = arma::dot(momentum,electric_field);

  // Set the momentum differentials.
  dxdt(4)                   = charge_to_rel_mass*arma::dot(momentum,electric_field);
  dxdt.subvec(5,7)          = charge*electric_field+charge_to_rel_mass*arma::cross(momentum,magnetic_field);

  // Model the effects of radiation reaction, if activated by the user.
  if (radiation_reaction != NoRR)
  {
    // Lorentz vector.
    arma::colvec lorentz = gamma*electric_field + arma::cross(momentum,magnetic_field);

    // Cheeky chi prefactor adds a dimensionful quantity to the mix.
    double chi_prefac        = constants::physics::hbar*unit_system.omega_0_SI/(mass*constants::physics::electron_mass*std::pow(constants::physics::c,2));
    chi_sq                   = std::pow(chi_prefac,2)*(std::pow(arma::norm(lorentz,2),2)-std::pow(pdotE,2));

    if (radiation_reaction == LandauLifshitz)
    {
      dxdt.subvec(4,7)      -= 2.0*constants::physics::alpha*std::pow(charge/mass,4)/(3.0*gamma)*chi_sq/chi_prefac*x.subvec(4,7);
    }

    else if (radiation_reaction == LandauLifshitzQuantumCorrection)
    {
      double quantum_factor  = std::pow(1+18.0*chi+69.0*chi_sq*73.0*std::pow(chi,3)+5.806*std::pow(chi_sq,2),-1.0/3.0);
      dxdt.subvec(4,7)      -= quantum_factor*2.0*constants::physics::alpha*std::pow(charge/mass,4)/(3.0*gamma)*chi_sq/chi_prefac*x.subvec(4,7);
    }
  }
}

template <class FieldModel>
inline
ParticleIonized<FieldModel>::ParticleIonized(const double                     my_charge,
                                             const double                     my_mass,
                                                   FieldModel              &  my_field_model,
                                                   MellotronUnits          &  my_units,
                                                   double                     my_field_threshold,
                                             const RadiationReactionModel     my_radiation_reaction)
: Particle<FieldModel>(my_charge,my_mass,my_field_model,my_units,my_radiation_reaction)
, field_threshold(my_field_threshold)
, ApplyLorentzForce(false)
{}

template <class FieldModel>
inline
void
ParticleIonized<FieldModel>::operator() (const arma::colvec::fixed<8>  &  x,
                                               arma::colvec::fixed<8>  &  dxdt,
                                         const double                     t)
{
  // We compute the field tensor.
  this->ComputeFieldTensor(x[0],x[1],x[2],x[3]);

  // We check whether the field has attained the threshold.
  if (!ApplyLorentzForce)
  {
    ApplyLorentzForce = arma::norm(this->electric_field,2) > this->field_threshold;
  }

  if (ApplyLorentzForce)
    this->ComputeLorentzForce(x, dxdt, t);

  else
    dxdt.fill(0.0);
}

} // namespace mellotron

#endif // PARTICLE_MEAT_HPP

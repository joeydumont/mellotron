#ifndef E_DIPOLES_HPP
#define E_DIPOLES_HPP

#include <boost/math/constants/constants.hpp>

namespace cst = boost::math::constants;

int interface_to_cubature_e_dipoles(unsigned int ndim, const double * x,    void *fdata,
                                    unsigned int fdim,       double * fval);

/*!
 *  \class  DipoleQuasiGaussian
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2016-10-18
 *  \brief  Implements the linearly polarized Salamin model.
 *
 * Evaluates the electromagnetic field components of an e-dipole pulse with
 * a quasi-Gaussian envelope.
 *
 * Ref: I. Gonoskov, Phys. Rev. A. 86, 053836 (2012).
 */
class DipoleQuasiGaussian
{
public:

  /// The model depends on the central frequency of the driving function omega,
  /// the pulse duration a and the laser energy energy.
  DipoleQuasiGaussian(double my_omega,double my_pulse_duration,double my_energy)
  : omega(my_omega)
  , lambda(2.0*cst::pi<double>()/omega)
  , a(my_pulse_duration)
  , energy(my_energy)
  {
    SetD0();
  }

  /// Computes the field components.
  std::array<double,6> ComputeFieldComponents(double t, double x, double y, double z)
  {
    // Prepare return value.
    std::array<double,6> field = {0.0,0.0,0.0,0.0,0.0,0.0};

    // Auxiliary variables.
    double R = std::sqrt(x*x+y*y+z*z);

    if (R < 1e-2*lambda)
    {
      field[2] = 4.0*d*g3(t)/3.0;
      return field;
    }

    // Derivatives of the driving function.
    double my_g0m = g0m(t,R);
    double my_g1p = g1p(t,R);
    double my_g2m = g2m(t,R);
    double my_g2p = g2p(t,R);

    // Auxiliary variables.
    double rm1 = 1.0/R;
    double rm2 = 1.0/(R*R);
    double rm3 = 1.0/(R*R*R);


    // Auxiliary electric fields.
    double E1 = rm1*my_g2m;
    double E2 = rm2*my_g1p + rm3*my_g0m;

    // Actual fields.
    double Ex = d*x*z*rm2*(E1+3.0*E2);
    double Ey = d*y*z*rm2*(E1+3.0*E2);
    double Ez = d*((-x*x-y*y)*rm2*E1+(3.0*z*z*rm2-1.0)*E2);

    // Auxiliary magnetic field.
    double B1 = my_g2p*rm1+my_g2m*rm2;

    // Actual fields.
    double Bx = -d*y*rm1*B1;
    double By = d*x*rm1*B1;
    double Bz = 0.0;


    field = {Ex,Ey,Ez,Bx,By,Bz};

    return field;
  }


    /// Compute the energy contained in the field.
  int ComputeNormalizationFactor()
  {
    const uint ndim = 3;
    const uint fdim = 1;
    double xmin[3] = {-5.0*lambda,-5.0*lambda,-5.0*lambda};
    double xmax[3] = { 5.0*lambda, 5.0*lambda, 5.0*lambda};
    double val[1], err[1];

    int error_flag = hcubature(fdim, interface_to_cubature_e_dipoles, this, ndim, xmin, xmax,
                               0,0,1.0e-6, ERROR_INDIVIDUAL, val, err);

    return error_flag;
  }

protected:
  void SetD0()
  {
    double sqr = std::sqrt(3.0*a*energy/(4.0*cst::pi<double>()))/(omega*omega);
    double xi  = a/omega;
    double f   = std::pow(1.0+6.0*xi*xi+3.0*std::pow(xi,4)*(1.0-std::exp(-0.5/(xi*xi))),-0.5)
                *std::pow(cst::pi<double>()/2.0, -0.25);

    d = sqr*f;
  }

  double g0_ret(double t, double R)
  {
    return std::exp(-a*a*std::pow(t-R,2))*std::sin(omega*(t-R));
  }

  double g0_adv(double t, double R)
  {
    return std::exp(-a*a*std::pow(t+R,2))*std::sin(omega*(t+R));
  }

  double g0p(double t, double R)
  {
    return g0_ret(t,R)+g0_adv(t,R);
  }

  double g0m(double t, double R)
  {
    return g0_ret(t,R)-g0_adv(t,R);
  }

  double g1p(double t, double R)
  {
    return (-2.0*a*a*(t+R)*std::sin(omega*(t+R))+std::cos(omega*(t+R))*omega)*std::exp(-std::pow(a*(t+R),2))
          +(-2.0*a*a*(t-R)*std::sin(omega*(t-R))+std::cos(omega*(t-R))*omega)*std::exp(-std::pow(a*(t-R),2));
  }

  double g1m(double t, double R)
  {
    return (2.0*a*a*(t+R)*std::sin(omega*(t+R))-std::cos(omega*(t+R))*omega)*std::exp(-std::pow(a*(t+R),2))
         +(-2.0*a*a*(t-R)*std::sin(omega*(t-R))+std::cos(omega*(t-R))*omega)*std::exp(-std::pow(a*(t-R),2));
  }

  double g2p(double t, double R)
  {
    return (-2.0*a*a*std::sin(omega*(t+R))+4.0*std::pow(a*a*(t+R),2)*std::sin(omega*(t+R))-4.0*a*a*(t+R)*std::cos(omega*(t+R))*omega-std::sin(omega*(t+R))*omega*omega)*std::exp(-std::pow(a*(t+R),2))
          +(-2.0*a*a*std::sin(omega*(t-R))+4.0*std::pow(a*a*(t-R),2)*std::sin(omega*(t-R))-4.0*a*a*(t-R)*std::cos(omega*(t-R))*omega-std::sin(omega*(t-R))*omega*omega)*std::exp(-std::pow(a*(t-R),2));
  }

  double g2m(double t, double R)
  {
    return (2.0*a*a*std::sin(omega*(t+R))-4.0*std::pow(a*a*(t+R),2)*std::sin(omega*(t+R))+4.0*a*a*(t+R)*std::cos(omega*(t+R))*omega+std::sin(omega*(t+R))*omega*omega)*std::exp(-std::pow(a*(t+R),2))
         +(-2.0*a*a*std::sin(omega*(t-R))+4.0*std::pow(a*a*(t-R),2)*std::sin(omega*(t-R))-4.0*a*a*(t-R)*std::cos(omega*(t-R))*omega-std::sin(omega*(t-R))*omega*omega)*std::exp(-std::pow(a*(t-R),2));
  }

  double g3(double t)
  {
    double exp_prefac = std::exp(-a*a*t*t);
    double cos_prefac = -6.0*a*a*omega+12*std::pow(a,4)*omega*t*t-std::pow(omega,3);
    double sin_prefac = 12.0*std::pow(a,4)*t-8.0*std::pow(a,6)*std::pow(t,3)+6.0*a*a*omega*omega*t;
    return exp_prefac*(cos_prefac*std::cos(omega*t)+sin_prefac*std::sin(omega*t));
  }

  double omega;
  double lambda;
  double a;
  double energy;
  double d;

};

/// Interface to Cubature.
int interface_to_cubature_e_dipoles(unsigned int ndim, const double * x,    void *fdata,
                                    unsigned int fdim,       double * fval)
{

  // We compute the electromagnetic field.
  auto obj = (DipoleQuasiGaussian * )fdata;
  auto field = obj->ComputeFieldComponents(0.0, x[0], x[1], x[2]);

  // We compute the electromagnetic energy.
  fval[0] = 0.0;
  for (uint i=0; i<6; i++)
  {
    fval[0] += field[i]*field[i];
  }

  fval[0] *= 0.5;

  return 0;
}

#endif // E_DIPOLES_HPP

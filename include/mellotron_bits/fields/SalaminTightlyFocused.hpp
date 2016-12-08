#ifndef SALAMIN_TIGHTLY_FOCUSED_HPP
#define SALAMIN_TIGHTLY_FOCUSED_HPP

#include <Cubature>
#include <boost/math/constants/constants.hpp>
#include <complex>

using namespace std::complex_literals;

namespace mellotron {
namespace cst = boost::math::constants;

/// Forward declaration of the interface to Cubature.
int interface_to_cubature_salamin(unsigned int ndim, const double * x,    void *fdata,
                                  unsigned int fdim,       double * fval);

/*!
 *  \class  SalaminTightlyFocusedLinear
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2016-10-18
 *  \brief  Implements the linearly polarized Salamin model.
 *
 * Evaluates the electromagnetic field components of the Salamin solution, linearly polarized
 * case. Some care must be taken near the origin, as the some of the auxiliary functions diverge.
 *
 * Ref: Y. I. Salamin, Phys. Rev. A. 92, 063818 (2015).
 */
class SalaminTightlyFocusedLinear
{
public:

  /// The model depends on the following parameters, the wavelength lambda,
  /// the waist waist and the axial length L.
  SalaminTightlyFocusedLinear(double my_lambda,double my_waist,double my_L,double my_energy)
  : lambda(my_lambda)
  , waist(my_waist)
  , L(my_L)
  , energy(my_energy)
  {
    norm_factor = 1.0;
    ComputeNormalizationFactor();
  }

  /// Computes the field components.
  std::array<double,6> ComputeFieldComponents(double t, double x, double y, double z)
  {
    double pi    = cst::pi<double>();
    // Define the beam coordinates.
    double r     = std::sqrt(x*x+y*y);
    double rho   = r / waist;
    double eta   = 0.5*(z+t);
    double zeta  = z - t;
    double zetap = pi*zeta/L;

    // Define some beam properties.
    double               k0    = 2.0*pi/lambda;
    double               zr    = pi*std::pow(waist,2)/lambda;

    // Define auxiliary variables.
    double               alpha  = eta/zr;
    std::complex<double> p      = std::complex<double>(1.0,alpha);
    std::complex<double> E0     = 2.0*k0*zr/(1i*(1.0+4.0*k0*zr));
    std::complex<double> Q1     = 4.0*1i*k0 - 2.0/zeta + 1i*(p-rho*rho)/(zr*p*p);
    std::complex<double> Q2     = 4.0*1i*k0 - 2.0/zeta - 1i*(2.0*p-rho*rho)/(zr*p*p);
    std::complex<double> Q3     = 4.0/(zeta*zeta) + (p-2.0*rho*rho)/(zr*zr*p*p*p) - std::pow((2.0*pi/L)/std::sin(zetap),2);
    std::complex<double> R      = 0.5*(Q1 + 2.0*pi/(L*std::tan(zetap)));
    std::complex<double> prefac = E0/k0*std::exp(2.0*1i*k0*zeta)/(zetap)*std::exp(-rho*rho/p)/p;

    // Compute the fields.
    double Ex = std::real(0.5*prefac
                          *(
                            (Q1 - (16.0*x*x/(p*p*std::pow(waist,4)))*(1.0/(2.0*R)-1i/(4.0*zr*p*R*R)))*std::sin(zetap)
                            +2.0*pi/L*std::cos(zetap)
                           )
                         );


    double Ey = std::real(-8.0*prefac*x*y/(p*p*std::pow(waist,4))*(1.0/(2.0*R)-1i/(4.0*zr*p*R*R))*std::sin(zetap));
    double Ez = std::real(2.0*prefac*x/(p*waist*waist)
                          *(
                            (Q2/(2.0*R)-Q3/(4.0*R*R))*std::sin(zetap)
                            +pi/(R*L)*std::cos(zetap)
                            )
                          );
    double Bx = 0.0;
    double By = std::real(prefac/2.0
                          *(
                            (4.0*1i*k0-2.0/zeta)*std::sin(zetap)
                            +2.0*pi/L*std::cos(zetap)
                            )
                          );
   double Bz = std::real(2.0*prefac*y/(p*waist*waist)*std::sin(zetap));

   // Determine limits by hand.
   if (std::abs(zetap) < 1.0e-45)
   {
    Ex = std::real(prefac*pi/L);
    Ey = 0.0;
    Ez = std::real(2.0*prefac*x/(p*waist*waist)*pi/(R*L));
    Bx = 0.0;
    By = std::real(prefac*pi/L);
    Bz = 0.0;
   }

   std::array<double,6> field = {Ex,Ey,Ez,Bx,By,Bz};
   for (uint i=0; i<6; i++)
   {
    field[i] /= norm_factor;
   }


   return field;
  }

  double lambda;
  double waist;
  double L;
  double energy;
  double norm_factor;

protected:
  /// Compute the energy contained in the field and rescale the components.
  int ComputeNormalizationFactor()
  {
    const uint ndim = 3;
    const uint fdim = 1;
    double xmin[3] = {-10.0*lambda,-10.0*lambda,-10.0*L};
    double xmax[3] = {10.0*lambda,10.0*lambda,10.0*L};
    double val[1], err[1];

    int error_flag = hcubature(fdim, interface_to_cubature_salamin, this, ndim, xmin, xmax,
                               0,0,1.0e-5, ERROR_INDIVIDUAL, val, err);

    norm_factor = std::sqrt(val[0]/energy);

    return error_flag;
  }

};

/// Interface to Cubature.
int interface_to_cubature_salamin(unsigned int ndim, const double * x,    void *fdata,
                                  unsigned int fdim,       double * fval)
{

  // We compute the electromagnetic field.
  auto obj = (SalaminTightlyFocusedLinear * )fdata;
  auto field = obj->ComputeFieldComponents(0.2*obj->lambda, x[0], x[1], x[2]);

  // We compute the electromagnetic energy.
  fval[0] = 0.0;
  for (uint i=0; i<6; i++)
  {
    fval[0] += field[i]*field[i];
  }

  fval[0] *= 0.5;

  return 0;
}

} // namespace mellotron

#endif // SALAMIN_TIGHTLY_FOCUSED_HPP

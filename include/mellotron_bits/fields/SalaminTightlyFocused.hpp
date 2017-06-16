#ifndef SALAMIN_TIGHTLY_FOCUSED_HPP
#define SALAMIN_TIGHTLY_FOCUSED_HPP

#include <Cubature>
#include <cuba.h>
#include <boost/math/constants/constants.hpp>
#include <complex>
#include <stdexcept>

using namespace std::complex_literals;

namespace mellotron {
namespace cst = boost::math::constants;

/// enum to choose the quadrature method.
enum SalaminQuadratureMethod {CubatureH, CubatureP, CubaCuhre};

/// Forward declaration of the interface to Cubature.
int interface_to_cubature_salamin(unsigned int ndim, const double * x,    void *fdata,
                                  unsigned int fdim,       double * fval);

/// Forward declaration of the interface to Cuba.
int interface_to_cuba_salamin(const int *ndim, const double x[],
	                            const int *fdim, double fval[], void *fdata);

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
  /// the waist waist and the axial length L. It is also possible to pass
  /// parameters for the integration routine.
  ///		@param[in] my_lambda		Central wavelength of the beam.
  ///   @param[in] my_waist			Waist size of the beam.
  ///   @param[in] my_L         Axial length of the beam.
  ///   @param[in] my_method    Integration method to use.
  ///   @param[in] my_xmin      Lower bounds of the integration. Array in units of (lambda, lambda, L).
  ///   @param[in] my_xmax      Upper bounds of the integration. Array in units of (lambda, lambda, L).
  ///   @param[in] my_abstol    Absolute tolerance.
  ///   @param[in] my_reltol    Relative tolerance.
  ///   @param[in] my_mineval   Minimum number of function evaluations.
  ///   @param[in] my_maxeval   Maximum number of function evaluations.
  SalaminTightlyFocusedLinear(double my_lambda,double my_waist,double my_L,double my_energy,
  	SalaminQuadratureMethod my_method = CubaCuhre, std::initializer_list<double> my_xmin = {-15,-15,-30}, std::initializer_list<double> my_xmax = {15,15,30},
  	double my_abstol = 1.0e-4, double my_reltol = 1.0e-4, int my_mineval = 1, int my_maxeval = 1e8)
  : lambda(my_lambda)
  , waist(my_waist)
  , L(my_L)
  , energy(my_energy)
  , method(my_method)
  , xmin(my_xmin)
  , xmax(my_xmax)
  , abstol(my_abstol)
  , reltol(my_reltol)
  , mineval(my_mineval)
  , maxeval(my_maxeval)
  {
  	xmin[0] *= lambda;
   	xmin[1] *= lambda;
   	xmin[2] *= L;
   	xmax[0] *= lambda;
   	xmax[1] *= lambda;
   	xmax[2] *= L;

    norm_factor = 1.0;
    ComputeNormalizationFactor();
  }

  /// This constructor takes an additional argument, my_norm_constant, which consists in the integral
  /// of 0.5*(E^2+B^2) for an energy of 1.0.
  SalaminTightlyFocusedLinear(double my_lambda,double my_waist,double my_L,double my_norm_constant,double my_energy)
  : lambda(my_lambda)
  , waist(my_waist)
  , L(my_L)
  , energy(my_energy)
  {
    norm_factor = my_norm_constant/std::sqrt(energy);
  }

  /// Computes the field components.
  std::array<double,6> ComputeFieldComponents(double t, double x, double y, double z) const
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
    std::complex<double> Q3     = 4.0/(zeta*zeta) + (p-2.0*rho*rho)/(zr*zr*p*p*p) - std::pow(2.0*pi/(L*std::sin(zetap)),2);
    std::complex<double> R      = 0.5*(Q1 + 2.0*pi/(L*std::tan(zetap)));
    std::complex<double> prefac = E0/k0*std::exp(2.0*1i*k0*zeta)/(zetap)*std::exp(-rho*rho/p)/p;

    // Compute the fields.
    double Ex = std::real(0.5*prefac
                          *(
                            (Q1 - (16.0*x*x/(p*p*std::pow(waist,4)))*(0.5/R-1i/(4.0*zr*p*R*R)))*std::sin(zetap)
                            +2.0*pi/L*std::cos(zetap)
                           )
                         );


    double Ey = std::real(-8.0*prefac*x*y/(p*p*std::pow(waist,4))*(0.5/R-1i/(4.0*zr*p*R*R))*std::sin(zetap));
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

  double lambda;           ///< Central wavelength of the pulse.
  double waist;            ///< Beam waist size (radius).
  double L;                ///< Axial length of the pulse, related to FWHM duration.
  double energy;           ///< Total energy contained in the pulse.
  double norm_factor;      ///< Normalization factor for the field to contain the proper energy.

  /// Compute the energy contained in the field and rescale the components.
  int ComputeNormalizationFactor()
  {
  	const int ndim = 3;
  	const int fdim = 1;
    int nregions, neval, error_flag;
    double val[1], err[1], prob[1];

    if (method == CubaCuhre)
    {
	    Cuhre(ndim,                           // Number of dimensions over which to integrate.
	    		  fdim,                           // Number of components of the vector integrand. We have a scalar integrand here.
	    		  interface_to_cuba_salamin,      // Function that evaluates the integrand.
	    		  this,                           // Pass the pointer to the object currently allocated. Allows Cuhre to access elements inside the class.
	    		  1,                              // Number of points to be given to the integrand for evaluation. Not sure how to make ComputeFieldComponents SIMD-compatible, so set to 1.
	    		  reltol,abstol,                  // Relative and absolute tolerance.
	    		  0,                              // Complicated system of flags passed to Cuhre and other integrators. Unused here.
	    		  mineval,maxeval,                // The minimum and maximum number of function evaluations.
	    		  0,                              // Specifies the order of the integration. In 3D, this results in degree-11 rule.
	    		  NULL,                           // Can be used to specify a checkpoint/restart file.
	    		  NULL,                           // Spin variable, used in the threading model. Best not to use it.
	    	    &nregions,&neval,&error_flag,		// Stores the number of regions that were needed, the number of function evaluations that were needed and the error flag.
	    	    &val[0],&err[0],&prob[0]);      // Stores the actual value of the integration, the estimated error and the xi^2 probability of error (not really applicable to Cuhre).
    }

    else if (method == CubatureH)
    {
    	error_flag = hcubature(fdim,                              // Number of components of the vector integrand. We have a scalar integrand here.
    		                     interface_to_cubature_salamin,     // Function that evaluates the integrand.
    		                     this,															// Passes the pointer to the object currently allocated. Allows Cuhre to access elements inside the class.
    		                     ndim,                              // Number of dimensions over which to integrate.
    		                     xmin.data(), xmax.data(),          // Integration boundaries.
    	                       maxeval,                           // Maximum number of evaluations.
    	                       abstol,reltol, ERROR_INDIVIDUAL,		// Absolute and relative tolerances, and the error norm (irrelevant as fdim=1).
    	                       val, err);                         // Store the value of the integrand and its error.
    }

    else if (method == CubatureP)
    {
    	error_flag = pcubature(fdim,                              // Number of components of the vector integrand. We have a scalar integrand here.
    		                     interface_to_cubature_salamin,     // Function that evaluates the integrand.
    		                     this,															// Passes the pointer to the object currently allocated. Allows Cuhre to access elements inside the class.
    		                     ndim,                              // Number of dimensions over which to integrate.
    		                     xmin.data(), xmax.data(),          // Integration boundaries.
    	                       maxeval,                           // Maximum number of evaluations.
    	                       abstol,reltol, ERROR_INDIVIDUAL,		// Absolute and relative tolerances, and the error norm (irrelevant as fdim=1).
    	                       val, err);                         // Store the value of the integrand and its error.
    }

    else
    {
    	throw std::runtime_error("Unknown integration method in SalaminTightlyFocusedLinear");
    }

    if (error_flag > 0)
    {
    	std::cerr << "There was an error in ComputeNormalizationFactor: " << error_flag << "\n";
    	std::cerr << "The result is " << val[0] << " ± " << err[0] << "." << "\n";
    }

    norm_factor = std::sqrt(val[0]/energy);

    return error_flag;
  }

	SalaminQuadratureMethod		method;
	std::vector<double>       xmin;
	std::vector<double>       xmax;
	double 										abstol;
	double 										reltol;
	int 											mineval;
	int 											maxeval;

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

/// Interface to Cuba
int interface_to_cuba_salamin(const int *ndim, const double x[],
	                            const int *fdim, double fval[], void *fdata)
{
	auto obj = (SalaminTightlyFocusedLinear* ) fdata;

  // We scale the integrand by applying the transformation
  // 	int_a^b f(x)dx --> int_0^1 f[a+(b-a)*y]*(b-a)dy
  // in each dimension.

  double xscaled[3];
  double jacobian = 1.0;
  for (int i=0; i<(*ndim); i++)
  {
  	double range = obj->xmax[i]-obj->xmin[i];
  	jacobian *= range;
  	xscaled[i] = obj->xmin[i]+range*x[i];
  }

	// We compute the electromagnetic field.
	auto field = obj->ComputeFieldComponents(0.2*obj->lambda, xscaled[0],xscaled[1], xscaled[2]);

	// We compute the electromagnetic energy.
	fval[0] = 0.0;
	for (uint i=0; i<6; i++)
	{
		fval[0] += field[i]*field[i];
	}

	fval[0] *= 0.5*jacobian;

	return 0;
}

} // namespace mellotron

#endif // SALAMIN_TIGHTLY_FOCUSED_HPP

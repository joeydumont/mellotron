/*! ------------------------------------------------------------------------- *
 * \author Denis Gagnon 										              *
 * \since 2016-11-29                                                          *
 *                                                                            *
 * Simulation program for election/ion trajectory computation using           *
 * Salamin's model. Reads parameters from the command line.		              *
 * --------------------------------------------------------------------------*/

#include <armadillo>
#include <mellotron>
#include <boost/program_options.hpp>
#include <boost/numeric/odeint.hpp>

namespace po = boost::program_options;
namespace odeint = boost::numeric::odeint;

int main(int argc, char* argv[])
{

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help", "produce help message")
    ("init_conds", po::value<std::vector<double> >()->multitoken(), "Initial position and momentum, respectively (6-vector")
    ("energy",     po::value<double>()->default_value(1.0),         "Pulse energy (joules)"                      )
    ("lam",        po::value<double>()->default_value(0.8),         "Wavelength in microns"                      )
    ("w0",         po::value<double>()->default_value(1.0),         "Beam waist in units of wavelength"          )
    ("L",          po::value<double>()->default_value(1.0),         "Axial length of beam in units of wavelength")
    ("mass",       po::value<double>()->default_value(1.0),         "Particle mass in units of electron mass"    )
    ("Q",          po::value<double>()->default_value(-1.0),        "Particle charge in units of electron charge")
    ("inttime",    po::value<double>()->default_value(6.28),        "Total integration time"                     )
    ("steps",      po::value<double>()->default_value(1e-01),       "Number of time steps"                       )
    ;


    // Parse command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    
    // Control the number of components in initial conditions vector
    std::vector<double> init_conds;
    if (!vm["init_conds"].empty() && (init_conds = vm["init_conds"].as<std::vector<double> >()).size() == 7)
    {
        // Good to go
    }
    else
    {
        std::cout 
        << "Initial conditions must be a 7 component vector... Exiting." 
        << std::endl;
        return 0;
    }
    
    // Assign other values to variables
    double energy  = vm["energy"].as<double>();
    double lam     = vm["lam"].as<double>();
    double w0      = vm["w0"].as<double>();
    double L       = vm["L"].as<double>();
    double mass    = vm["mass"].as<double>();
    double Q       = vm["Q"].as<double>();
    double tmin    = vm["tmin"].as<double>();
    double inttime = vm["inttime"].as<double>();
    double dt      = vm["dt"].as<double>();
    
    // Create field object
    SalaminTightlyFocusedLinear                  field(lam,w0,L,energy);     
    Particle<SalaminTightlyFocusedLinear>             ion(Q,mass,field);
    ParticleObserver<SalaminTightlyFocusedLinear>          ion_obs(ion);
    
    // Define the initial conditions.
    arma::colvec::fixed<8> x = arma::zeros<arma::colvec>(8);
    x[0] = init_conds[0]; x[1] = init_conds[1]; x[2] = init_conds[2]; x[3] = init_conds[3];
    //x[0] = init_conds[1];

    return 0;

}

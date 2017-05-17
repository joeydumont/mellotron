/*! ------------------------------------------------------------------------- *
 * \author Denis Gagnon 										              *
 * \since 2016-11-29                                                          *
 *                                                                            *
 * Simulation program for election/ion trajectory computation using           *
 * Salamin's model. Reads parameters from the command line.		              *
 * --------------------------------------------------------------------------*/

#include <armadillo>
#include <cmath>
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
    ("init_conds", po::value<std::vector<double> >()->multitoken(), "Initial position and momentum, electronic units (6-vector)")
    ("energy",     po::value<double>()->default_value(10.0)         "Pulse energy, in joules"                               )
    ("lam",        po::value<double>()->default_value(8.0e-07),     "Wavelength in meters"                                  )
    ("w0",         po::value<double>()->default_value(1.0),         "Beam waist in units of lambda"                         )
    ("L",          po::value<double>()->default_value(1.0),         "Axial length of beam in units of lambda"               )
    ("mass",       po::value<double>()->default_value(1.0),         "Particle mass in units of electron mass"               )
    ("Q",          po::value<double>()->default_value(-1.0),        "Particle charge in units of electron charge"           )
    ("t_init",     po::value<double>()->default_value(-5.0e-12),    "Initial time in simulation in seconds"           )
    ("dt",         po::value<double>()->default_value(1e-13),       "Duration of a time step in seconds"              )
    ("nsteps",     po::value<int>()->default_value(100),            "Number of time steps"                            )
    ;

    // Parse command line and store in variable map
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // Control the number of components in initial conditions vector
    std::vector<double> init_conds;
    if (!vm["init_conds"].empty() && (init_conds = vm["init_conds"].as<std::vector<double> >()).size() == 6)
    {
        // Good to go
        init_conds = vm["init_conds"].as<std::vector<double> >();
    }
    else
    {
        std::cout
                << "Initial conditions must be a 6 component vector... Exiting."
                << std::endl;
        return 0;
    }

    // Parse lambda from command line
    double lam     = vm["lam"].as<double>();
    
    // Instantiate electron units object
    mellotron::MellotronUnits electron_units
    (2.0*mellotron::constants::math::pi*mellotron::constants::physics::c/lam);
    
    // Convert everything to electronic units (get everything from vm, except lam)
    lam /= electron_units.UNIT_LENGTH;
    
    double energy  = vm["energy"].as<double>() / electron_units.UNIT_ENERGY ;
    double w0      = vm["w0"].as<double>() * lam;
    double L       = vm["L"].as<double>() * lam;
    double mass    = vm["mass"].as<double>();
    double Q       = vm["Q"].as<double>();
    double t_init  = vm["t_init"].as<double>() / electron_units.UNIT_TIME ;
    double dt      = vm["dt"].as<double>() / electron_units.UNIT_TIME ;
    double nsteps  = vm["nsteps"].as<int>();

    // Create field object
    mellotron::SalaminTightlyFocusedLinear                              field(lam,w0,L,energy);
    mellotron::Particle<mellotron::SalaminTightlyFocusedLinear>         particle(Q,mass,field);
    mellotron::ParticleObserver<mellotron::SalaminTightlyFocusedLinear> particle_obs(particle);

    // Define the initial conditions.
    arma::colvec::fixed<8> x = arma::zeros<arma::colvec>(8);
    double x_init                      = init_conds[0]; // Initial position
    double y_init                      = init_conds[1];
    double z_init                      = init_conds[2];

    double px_init = init_conds[3]; // Initial momentum
    double py_init = init_conds[4];
    double pz_init = init_conds[5];


    // Times at which we output the data.
    arma::colvec times = arma::linspace<arma::colvec>(t_init,t_init+nsteps*dt,nsteps); // Time vector

    // Set the initial conditions.
    particle.SetInitConditions(x,x_init,y_init,z_init,px_init,py_init,pz_init,times[0]);

    // Perform integration
    size_t steps = odeint::integrate_times(
                       odeint::make_dense_output(1.0e-6,1.0e-6, odeint::runge_kutta_dopri5<arma::colvec::fixed<8> >() ),
                       std::ref(particle),
                       x,
                       boost::begin(times),
                       boost::end(times),
                       0.1,
                       std::ref(particle_obs)
                   );

    std::cout << steps << std::endl;

    for (uint i=0; i<8; i++)
    {
        std::cout << x[i] << std::endl;
    }

    particle_obs.OutputData();

    return 0;

}

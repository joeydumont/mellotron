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
#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/numeric/odeint.hpp>

#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

using boost::math::relative_difference;
namespace po = boost::program_options;
namespace odeint = boost::numeric::odeint;

struct IntegrationSalaminConfig
{
    double lam_;
    double w0_;
    double L_;
    double energy_;
    double mass_;
    double Q_;
    std::string rad_react_;       // Radiation reaction model.
    double t_init_;
    double dt_;
    int nsteps_;
    void read(std::ifstream& file, IntegrationSalaminConfig*& config);
    std::string envelope_;
};

void IntegrationSalaminConfig::read(std::ifstream& file, IntegrationSalaminConfig*& config)
{
    using boost::property_tree::ptree;
    ptree pt;
    read_xml(file, pt);
    bool configIsEmpty = true;
    bool hasFoundParticle = false;
    bool hasFoundIntegSala = false;
    BOOST_FOREACH(ptree::value_type const& v, pt.get_child("config"))
    {
        if(v.first == "integration_salamin")
        {
            if(configIsEmpty)
            {
                config = new IntegrationSalaminConfig();
                configIsEmpty = false;
            }
            hasFoundIntegSala = true;
            config->lam_ = v.second.get<double>("lambda");
            config->w0_ = v.second.get<double>("w0");
            config->L_ = v.second.get<double>("L");
            config->energy_ = v.second.get<double>("energy");
            config->t_init_ = v.second.get<double>("t_init");
            config->dt_ = v.second.get<double>("dt");
            config->nsteps_ = v.second.get<int>("nsteps");
            config->envelope_ = v.second.get("envelope", "NoEnvelope");
        }
        if(v.first == "particle")
        {
            if(configIsEmpty)
            {
                config = new IntegrationSalaminConfig();
                configIsEmpty = false;
            }
            hasFoundParticle = true;
            config->mass_ = v.second.get<double>("mass");
            config->Q_ = v.second.get<double>("Q");
            config->rad_react_ = v.second.get("radiation_reaction", "NoRR");
        }
    }

    if(config == nullptr || !hasFoundParticle || !hasFoundIntegSala)
    {
        throw std::runtime_error("Missing integration_salamin or particle config.");
    }
}

int main(int argc, char* argv[])
{
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help",                                                                "produce help message"                                   )
    ("init_conds",         po::value<std::vector<double> >()->multitoken(), "Initial position and momentum, electronic units (6-vector)")
    ;

    // Parse command line and store in variable map
    po::variables_map vm;
    po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);

    // Open config file
    std::ifstream conf_file;
    conf_file.open("configSalamin.xml");
    if(!conf_file.is_open())
    {
        std::cout
                << "the config file must be in the same directory... Exiting."
                << std::endl;
        return 0;
    }

    // Read config file
    IntegrationSalaminConfig* config = nullptr;
    config->IntegrationSalaminConfig::read(conf_file, config);

    // Control the number of components in initial conditions vector
    bool IsConstantInitialTime = true;
    std::vector<double> init_conds;
    if (!vm["init_conds"].empty() && (init_conds = vm["init_conds"].as<std::vector<double> >()).size() == 6)
    {
        // Good to go
        init_conds = vm["init_conds"].as<std::vector<double> >();
    }
    else if (!vm["init_conds"].empty() && (init_conds = vm["init_conds"].as<std::vector<double> >()).size() == 7)
    {
        // Varying initial time. Must read inital time from CL, not from XML file.
        IsConstantInitialTime = false;
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
    double lam     = config->lam_;

    // Instantiate electron units object
    mellotron::MellotronUnits electron_units
    (2.0*mellotron::constants::math::pi*mellotron::constants::physics::c/lam);

    // Convert everything to electronic units (get everything from vm, except lam)
    lam /= electron_units.UNIT_LENGTH;

    double energy  = config->energy_ / electron_units.UNIT_ENERGY ;
    double w0      = config->w0_ * lam;
    double L       = config->L_ * lam;
    double mass    = config->mass_;
    double Q       = config->Q_;
    config->t_init_= config->t_init_ / electron_units.UNIT_TIME ;
    double dt      = config->dt_ / electron_units.UNIT_TIME ;
    int nsteps     = config->nsteps_;

    // Initial time
    // Behaviour depends on where we read the initial time.
    double t_init;
    if (IsConstantInitialTime)
        t_init = config->t_init_;
    else
        t_init = init_conds[6];

    // We verify that the normalization constant was calculated for the same
    // (lambda,w0,L) tuple.
    std::ifstream norm_constant_file;
    norm_constant_file.open("normalization_constant.txt");
    if(!norm_constant_file.is_open())
    {
        std::cout
                << "normalization_constant.txt must be in the same directory... Exiting."
                << std::endl;
        return 0;
    }

    std::string line;
    std::getline(norm_constant_file, line);std::getline(norm_constant_file, line);
    std::istringstream iss(line);
    double lambda_file,w0_file,L_file,norm_constant;
    iss >> lambda_file >> w0_file >> L_file >> norm_constant;

    if (
            relative_difference(config->lam_, lambda_file) > 1.0e-5
        ||  relative_difference(config->w0_,  w0_file)     > 1.0e-5
        ||  relative_difference(config->L_,   L_file)      > 1.0e-5)
    {
        throw std::runtime_error("Wrong value of the normalization constant.");
    }

    // Determine the radiation reaction model.
    mellotron::RadiationReactionModel rr_model;
    if (config->rad_react_ == std::string("NoRR"))
        rr_model = mellotron::NoRR;

    else if (config->rad_react_ == std::string("LandauLifshitz"))
        rr_model = mellotron::LandauLifshitz;

    else if (config->rad_react_ == std::string("LandauLifshitzQuantumCorrection"))
        rr_model = mellotron::LandauLifshitzQuantumCorrection;

    else
        rr_model = mellotron::NoRR;

    // Create envelope
    mellotron::Envelope* envelope;
    if (config->envelope_ == std::string("NoEnvelope"))
        envelope = new mellotron::NoEnvelope();

    else if (config->envelope_ == std::string("EnvelopeHann"))
        envelope = new mellotron::EnvelopeHann(t_init,config->nsteps_*dt);

    else
        envelope = new mellotron::NoEnvelope();

    // Create field object
    mellotron::SalaminTightlyFocusedLinear                              field(lam,w0,L,norm_constant,energy);
    mellotron::Particle<mellotron::SalaminTightlyFocusedLinear>         particle(Q,mass,field,electron_units,*envelope,rr_model);
    mellotron::ParticleObserver<mellotron::SalaminTightlyFocusedLinear> particle_obs(particle,nsteps);

    // Define the initial conditions.
    arma::colvec::fixed<8> x = arma::zeros<arma::colvec>(8);
    double x_init                      = init_conds[0]; // Initial position
    double y_init                      = init_conds[1];
    double z_init                      = init_conds[2];

    double px_init = init_conds[3]; // Initial momentum
    double py_init = init_conds[4];
    double pz_init = init_conds[5];


    // Times at which we output the data.
    arma::colvec times = arma::linspace<arma::colvec>(t_init,t_init+config->nsteps_*dt,config->nsteps_); // Time vector

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

    particle_obs.OutputData();

    delete envelope;
    delete config;
    return 0;

}

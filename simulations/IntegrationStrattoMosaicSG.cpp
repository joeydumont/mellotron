/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont      <joey.dumont@gmail.com>                           *
 * \since 2017-08-07                                                          *
 *                                                                            *
 * Simulation program using the strattocalculator via the                     *
 * StrattoCalculatorWrapper.hpp (mosaic fields).                              *
 * --------------------------------------------------------------------------*/

#include <cmath>
#include <meshpi>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/numeric/odeint.hpp>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "StrattoCalculatorWrapper.hpp"

using namespace MeshPI;
using namespace StrattoCalculator;
namespace po = boost::program_options;
namespace odeint = boost::numeric::odeint;

/// Parse an array in ini_parser.
/// http://stackoverflow.com/questions/4986052/boost-property-tree-working-with-simple-arrays-or-containers
template <typename T>
std::vector<T> to_array(const std::string &s)
{
  std::vector<T> result;
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss,item, ','))
    result.push_back(boost::lexical_cast<T>(item));
  return result;
}

// Structure in which we store all of the data we read in the config file.
struct StrattoMosaicConfig
{
    double r_min_;                // Minimum aperture
    double r_max_;                // Maximum aperture.
    double th_min_;               // Minimum angle
    double th_max_;               // Maximum angle (2*pi for a complete parabola)
    double focal_length_;         // Focal length of the parabola.
    bool boundary_min_;           // Flag to enable evaluation of boundary integral at r_min.
    bool boundary_max_;           // Flag to enable evaluation of boundary integral at r_max.
    unsigned int intervals_r_;    // Number of intervals in r coord
    unsigned int intervals_th_;   // Number of intervals in th coord
    unsigned int num_points_r_;   // Number of points per interval in r coord
    unsigned int num_points_th_;  // Number of points per interval in th coord
    double energy_;               // Energy in the incident beam
    double lambda_c_;             // Central wavelength of the spectrum
    double lambda_min_;           // Minimum wavelength of the spectrum
    double lambda_max_;           // Maximum wavelength of the spectrum
    double delta_lambda_;         // Full-width half-maximum of the spectrum
    int num_components_;          // Number of spectral components to consider
    int gaussian_order_;          // Order of the super-gaussian spectrum
    bool   hasChirp;              // Determines if the config file has chirp info.
    double omega_c_chirp_;         // Central frequency where we expand the spectral phase.
    std::vector<double> chirp_taylor_coefficients_;// Taylor coefficients of the phase (to model chirp).
    double beam_width_;           // 1/e radius of the field
    std::vector<double> lg_coeffs_;// Coefficients of the Laguerre-Gauss expansion.
    int    number_tesserae_;      // Number of segments in the mosaic.
    double mass_;                 // Particle mass
    double Q_;                    // Particle charge
    double t_init_;               // Initial time in simulation
    double dt_;                   // Duration of a time step
    int nsteps_;                  // Number of time steps
    void read(std::ifstream& file, StrattoMosaicConfig*& config);
};

// Function used to read the data read from the config file.
void StrattoMosaicConfig::read(std::ifstream& file, StrattoMosaicConfig*& config)
{
    using boost::property_tree::ptree;
    ptree pt;
    read_xml(file, pt);
    bool configIsEmpty = true;
    bool hasFoundParabola = false;
    bool hasFoundSpectrum = false;
    bool hasFoundModel = false;
    bool hasFoundParticle = false;
    bool hasFoundIntegration = false;
    BOOST_FOREACH(ptree::value_type const& v, pt.get_child("config"))
    {
        if(v.first == "parabola")
        {
            if(configIsEmpty)
            {
                config = new StrattoMosaicConfig();
                configIsEmpty = false;
            }
            hasFoundParabola = true;
            config->r_min_ = v.second.get<double>("r_min");
            config->r_max_ = v.second.get<double>("r_max");
            config->th_min_ = v.second.get<double>("th_min");
            config->th_max_ = v.second.get<double>("th_max");
            config->focal_length_ = v.second.get<double>("focal_length");
            config->boundary_min_ = v.second.get<bool>("boundary_min");
            config->boundary_max_ = v.second.get<bool>("boundary_max");
            config->intervals_r_ = v.second.get<int>("intervals_r");
            config->intervals_th_ = v.second.get<int>("intervals_th");
            config->num_points_r_ = v.second.get<int>("num_points_r");
            config->num_points_th_ = v.second.get<int>("num_points_th");
        }
        if(v.first == "spectrum")
        {
            if(configIsEmpty)
            {
                config = new StrattoMosaicConfig();
                configIsEmpty = false;
            }
            hasFoundSpectrum = true;
            config->energy_ = v.second.get<double>("energy");
            config->lambda_c_ = v.second.get<double>("lambda_c");
            config->lambda_min_ = v.second.get<double>("lambda_min");
            config->lambda_max_ = v.second.get<double>("lambda_max");
            config->delta_lambda_ = v.second.get<double>("delta_lambda");
            config->num_components_ = v.second.get<int>("num_components");
            config->gaussian_order_ = v.second.get<int>("gaussian_order");

            // Check that the config file contains chirp information.
            if (boost::optional<double> omega_c_chirp_opt = v.second.get_optional<double>("omega_c_chirp"))
            {
                config->hasChirp = true;
                config->omega_c_chirp_ = *omega_c_chirp_opt;
                config->chirp_taylor_coefficients_ = to_array<double>(v.second.get<std::string>("chirp_taylor_coefficients"));
            }

            else
                config->hasChirp = false;
        }
        if(v.first == "model")
        {
            if(configIsEmpty)
            {
                config = new StrattoMosaicConfig();
                configIsEmpty = false;
            }
            hasFoundModel = true;
            config->beam_width_ = v.second.get<double>("beam_width");
            config->lg_coeffs_  = to_array<double>(v.second.get<std::string>("lg_coeffs"));
            config->number_tesserae_ = v.second.get<int>("number_tesserae");
        }
        if(v.first == "particle")
        {
            if(configIsEmpty)
            {
                config = new StrattoMosaicConfig();
                configIsEmpty = false;
            }
            hasFoundParticle = true;
            config->mass_ = v.second.get<double>("mass");
            config->Q_ = v.second.get<double>("Q");
        }
        if(v.first == "integration")
        {
            if(configIsEmpty)
            {
                config = new StrattoMosaicConfig();
                configIsEmpty = false;
            }
            hasFoundIntegration = true;
            config->t_init_ = v.second.get<double>("t_init");
            config->dt_ = v.second.get<double>("dt");
            config->nsteps_ = v.second.get<int>("nsteps");
        }
    }

    if(config == nullptr || !hasFoundParabola || !hasFoundSpectrum || !hasFoundModel || !hasFoundParticle || !hasFoundIntegration)
    {
        throw std::runtime_error("Missing parabola, spectrum or model config.");
    }
}

double to_rad(double angle_in_deg)
{
  return angle_in_deg*GlobalConstant::PI/180.0;
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

    // Open config file.
    std::ifstream conf_file;
    conf_file.open("configStrattoMosaicSG.xml");
    if(!conf_file.is_open())
    {
        std::cout
                << "the config file must be in the same directory... Exiting."
                << std::endl;
        return 0;
    }

    // Read config file.
    StrattoMosaicConfig* config = nullptr;
    config->StrattoMosaicConfig::read(conf_file, config);

    // Instantiate electron units object.
    mellotron::MellotronUnits electron_units
    (2.0*mellotron::constants::math::pi*mellotron::constants::physics::c/config->lambda_c_);

    // Change the parabola parameters to appropriate units.
    config->r_min_        = config->r_min_/electron_units.UNIT_LENGTH;
    config->r_max_        = config->r_max_/electron_units.UNIT_LENGTH;
    config->th_min_       = to_rad(config->th_min_);
    config->th_max_       = to_rad(config->th_max_);
    config->focal_length_ = config->focal_length_/electron_units.UNIT_LENGTH;

    // Create the mesh of the parabola.
    std::array<double,2> min_parabola                      = {config->r_min_,config->th_min_};
    std::array<double,2> max_parabola                      = {config->r_max_,config->th_max_};
    std::array<unsigned int, 2> intervals_per_dim_parabola = {config->intervals_r_,config->intervals_th_};
    CoordinateSystemModel<2,Cylindrical> * coord_sys       = new CoordinateSystemModel<2,Cylindrical>(x12,space,qed_units);
    DomainSerial<2>* domain_parabola                       = new DomainSerial<2>(
        min_parabola,
        max_parabola,
        intervals_per_dim_parabola,
        coord_sys
    );

    std::array<unsigned int, 2> num_points_parabola        = {config->num_points_r_,config->num_points_th_};
    std::array<bool,2> boundary                            = {config->boundary_max_,config->boundary_min_};
    SurfaceMesh<2> *mesh_parabola                          = new SurfaceMesh<2>(
        domain_parabola,
        num_points_parabola,
        gausslegendre_mesh,
        boundary
    );
    EmittingSurfaceAxialSym<2> *surface                    = new EmittingSurfaceParabola<2>(
        mesh_parabola,
        config->focal_length_
        );
        surface->ComputeSurface();
        surface->ComputeNormal();

    // Change spectrum parameters to appropriate units.
    config->energy_       = config->energy_/electron_units.UNIT_ENERGY;
    config->lambda_c_     = config->lambda_c_/electron_units.UNIT_LENGTH;
    config->lambda_min_   = config->lambda_min_/electron_units.UNIT_LENGTH;
    config->lambda_max_   = config->lambda_max_/electron_units.UNIT_LENGTH;
    config->delta_lambda_ = config->delta_lambda_/electron_units.UNIT_LENGTH;

    // Create the spectrum.
    std::vector<double> beam_frequencies, beam_freq_widths, beam_freq_ratios;
    beam_frequencies.push_back(config->lambda_c_);
    beam_freq_widths.push_back(config->delta_lambda_);
    beam_freq_ratios.push_back(1.0);
    std::vector<double> phase(config->num_components_,0.0);
    Spectrum *spectrum_incident = new SpectrumGaussianWavelength(
        config->num_components_,
        config->lambda_min_,
        config->lambda_max_,
        phase,
        beam_frequencies,
        beam_freq_widths,
        beam_freq_ratios,
        config->energy_,
        config->gaussian_order_
    );
    for (unsigned int i=0; i<config->num_components_;i++)
    {
        spectrum_incident->SetPhase(i, -2.0*spectrum_incident->GetOmega(i)*config->focal_length_);
    }

    if (config->hasChirp)
    {
        double omega_c_chirp = config->omega_c_chirp_ * electron_units.UNIT_TIME;
        for (int i=0; i<config->chirp_taylor_coefficients_.size(); i++)
        {
            config->chirp_taylor_coefficients_[i] /= std::pow(electron_units.UNIT_TIME,i);
        }

        int factorial = 1;
        for (int i=0; i<config->num_components_; i++)
        {
            spectrum_incident->SetPhase(i, spectrum_incident->GetPhase(i)
                                            + config->chirp_taylor_coefficients_[i]*std::pow(spectrum_incident->GetOmega(i)-omega_c_chirp,i)/factorial);
            factorial *= (i+1);
        }
    }

    // Create the beam model.
    config->beam_width_ = config->beam_width_/electron_units.UNIT_LENGTH;
    TEM00ModeParaxial *beam = new TEM00ModeParaxial(spectrum_incident,config->beam_width_,config->lg_coeffs_,true);
    MosaicBeams<2>  *  mosaic = new MosaicBeams<2>(mesh_parabola,beam,config->number_tesserae_);

    // Evaluate the beam on the mirror.
    SurfaceEMFieldManyStorage<SurfaceEMFieldGeneral,2,1> *incident_field = new SurfaceEMFieldManyStorage<SurfaceEMFieldGeneral,2,1>(
        mesh_parabola,
        spectrum_incident,
        surface,
        mosaic
    );

    // Declare the integrator and convert it to the form the MELLOTRON expects.
    auto integrator
        = new StrattonChuIntegrator::MeshlessAxialSurfGenField(surface,incident_field);
    auto meshless_time_field
        = new TemporalEMFieldMeshless<StrattonChuIntegrator::MeshlessAxialSurfGenField>(integrator);
    auto wrapper
        = new StrattoCalculatorWrapper<StrattonChuIntegrator::MeshlessAxialSurfGenField>(*meshless_time_field);

    // Change the integration parameters to appropriate units.
    config->t_init_ = config->t_init_ / electron_units.UNIT_TIME;
    config->dt_ = config->dt_ / electron_units.UNIT_TIME;

    // MELLOTRON
    mellotron::Particle<StrattoCalculatorWrapper<StrattonChuIntegrator::MeshlessAxialSurfGenField>> particle(config->Q_,config->mass_,*wrapper,electron_units);
    mellotron::ParticleObserver<StrattoCalculatorWrapper<StrattonChuIntegrator::MeshlessAxialSurfGenField>> particle_obs(particle,config->nsteps_);

    // Define the initial conditions.
    arma::colvec::fixed<8> x = arma::zeros<arma::colvec>(8);
    double x_init                      = init_conds[0]; // Initial position
    double y_init                      = init_conds[1];
    double z_init                      = init_conds[2];

    double px_init = init_conds[3]; // Initial momentum
    double py_init = init_conds[4];
    double pz_init = init_conds[5];

    // Times at which we output the data.
    // Behaviour depends on where we read the initial time.
    double t_init;
    if (IsConstantInitialTime)
        t_init = config->t_init_;
    else
        t_init = init_conds[6];

    arma::colvec times = arma::linspace<arma::colvec>(t_init,t_init+config->nsteps_*config->dt_,config->nsteps_); // Time vector

    // Set the initial conditions.
    particle.SetInitConditions(x,x_init,y_init,z_init,px_init,py_init,pz_init,times[0]);

    // Perform integration.
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


    delete config;
    return 0;

}
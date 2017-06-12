/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont                                                        *
 * \since 2017-05-17                                                          *
 *                                                                            *
 * Computation of the normalization constant needed for the evaluation of the *
 * Salamin linearly polarized fields for a given set of parameters.           *
 * --------------------------------------------------------------------------*/

#include <armadillo>
#include <mellotron>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>


struct ComputeNormConstConfig
{
    double lam_;
    double w0_;
    double L_;
    void read(std::ifstream& file, ComputeNormConstConfig*& config);
};

void ComputeNormConstConfig::read(std::ifstream& file, ComputeNormConstConfig*& config)
{
    using boost::property_tree::ptree;
    ptree pt;
    read_xml(file, pt);
    BOOST_FOREACH(ptree::value_type const& v, pt.get_child("config"))
    {
        if(v.first == "normalization_constant")
        {
            config = new ComputeNormConstConfig();
            config->lam_ = v.second.get<double>("lambda");
            config->w0_ = v.second.get<double>("w0");
            config->L_ = v.second.get<double>("L");
        }
    }
    
    if(config == nullptr)
    {
        throw std::runtime_error("Missing normalization_constant config.");
    }
}


int main(int argc, char* argv[])
{
  // Open config file
  std::ifstream conf_file;
  conf_file.open("config.xml");
  if(!conf_file.is_open())
  {
      std::cout
          << "the config file must be in the same directory... Exiting."
          << std::endl;
      return 0;
  }

  // Read config file
  ComputeNormConstConfig* config = nullptr;
  config->ComputeNormConstConfig::read(conf_file, config);

  // Parse lambda to instantiate the proper MellotronUnits object.
  double lambda = config->lam_;

  // Compute the proper electronic units.
  mellotron::MellotronUnits electron_units(2.0*mellotron::constants::math::pi*mellotron::constants::physics::c/lambda);

  // Convert everything to electronic units.
         lambda /= electron_units.UNIT_LENGTH;

  double w0     = config->w0_ * lambda;
  double L      = config->L_  * lambda;

  // Create the field object.
  mellotron::SalaminTightlyFocusedLinear field(lambda,w0,L,1.0 / electron_units.UNIT_ENERGY );

  // Output the normalization constant in a file.
  std::ofstream outfile;
  outfile.open("normalization_constant.txt", std::ios::trunc);
  outfile << "# lambda \t"
          << "w0 \t\t"
          << "L \t\t"
          << "norm. constant \n";
  outfile << std::scientific << config->lam_
          << "\t"
          << std::scientific << config->w0_
          << "\t"
          << std::scientific << config->L_
          << "\t"
          << std::scientific << field.norm_factor
          << "\n";
  outfile.close();
  delete config;
  return 0;
}
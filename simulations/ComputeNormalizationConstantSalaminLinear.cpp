/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont                                                        *
 * \since 2017-05-17                                                          *
 *                                                                            *
 * Computation of the normalization constant needed for the evaluation of the *
 * Salamin linearly polarized fields for a given set of parameters.           *
 * --------------------------------------------------------------------------*/

#include <armadillo>
#include <mellotron>
#include <boost/program_options.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  // Declare the supported options.
  po::options_description desc("Mandatory arguments");
  desc.add_options()
  ("help",                                                    "produce help message")
  ("lambda",  po::value<double>()->default_value(800.0e-9),   "Wavelength of the beam in metres."  )
  ("w0",      po::value<double>()->default_value(0.8),        "Beam waist as a multiple of lambda.")
  ("L",       po::value<double>()->default_value(0.7),        "Axial length as a multiple of lambda.")
  ("outfile", po::value< std::string >()->required(),         "Name of the file to which we output the normalization constant")
  ;

  // Parse command line and store in variable map.
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  // Parse lambda to instantiate the proper MellotronUnits object.
  double lambda = vm["lambda"].as<double>();

  // Compute the proper electronic units.
  mellotron::MellotronUnits electron_units(2.0*mellotron::constants::math::pi*mellotron::constants::physics::c/lambda);

  // Convert everything to electronic units.
         lambda /= electron_units.UNIT_LENGTH;

  double w0     = vm["w0"].as<double>() * lambda;
  double L      = vm["L"].as<double>()  * lambda;

  // Create the field object.
  mellotron::SalaminTightlyFocusedLinear field(lambda,w0,L,1.0);

  // Output the normalization constant in a file.
  std::ofstream outfile;
  outfile.open((vm["outfile"].as<std::string>()).c_str(), std::ios::trunc);
  outfile << "# lambda \t"
          << "w0 \t\t"
          << "L \t\t"
          << "norm. constant \n";
  outfile << std::scientific << vm["lambda"].as<double>()
          << "\t"
          << std::scientific << vm["w0"].as<double>()
          << "\t"
          << std::scientific << vm["L"].as<double>()
          << "\t"
          << std::scientific << field.norm_factor
          << "\n";
  outfile.close();

  return 0;
}
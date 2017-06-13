/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont                  <joey.dumont@gmail.com>               *
 * \since 2016-10-28                                                          *
 *                                                                            *
 * Unit test of fields/eDipoles.cpp in MELLOTRON. We simply evaluate          *
 * and plot the field in a plane and compare with the original article.       *
 * --------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <armadillo>
#include <mellotron>
#include <boost/math/constants/constants.hpp>

using namespace mellotron;
class DipoleQuasiGaussianTest : public testing::Test
{
public:
  DipoleQuasiGaussianTest()
  : lambda(800e-9)
  , omega(2.0*constants::math::pi*constants::physics::c/lambda)
  , electronic_units(omega)
  {
    lambda         /= electronic_units.UNIT_LENGTH;
    omega          *= electronic_units.UNIT_TIME;
    pulse_duration = 0.1*omega;
    xmax           = lambda;
    energy         = 15.0/electronic_units.UNIT_ENERGY;
    field = new DipoleQuasiGaussian(omega,pulse_duration,energy);

    std::cout << "Energy [in Mellotron Units]: " << energy << std::endl;
  }

protected:

  virtual void SetUp()
  {}

  double lambda;
  double omega;
  MellotronUnits electronic_units;

  double pulse_duration;
  double xmax;
  double energy;

  DipoleQuasiGaussian  *  field;
};

TEST_F(DipoleQuasiGaussianTest, q_gauss)
{
  // Define the mesh of the plot.
  uint size_plot = 100;
  auto x_field = arma::linspace<arma::colvec>(-xmax,xmax, size_plot);
  auto y_field = arma::linspace<arma::colvec>(-xmax,xmax, size_plot);
  arma::mat Ex(size_plot,size_plot);
  arma::mat Ey(size_plot,size_plot);
  arma::mat Ez(size_plot,size_plot);
  arma::mat Bx(size_plot,size_plot);
  arma::mat By(size_plot,size_plot);
  arma::mat Bz(size_plot,size_plot);

  // Compute the field values.
  for (uint i=0; i<size_plot; i++)
  {
    for (uint j=0; j<size_plot; j++)
    {
      auto field_value  = field->ComputeFieldComponents(0.0*lambda,x_field[i],y_field[j],0.0);
      Ex(i,j) = field_value[0];
      Ey(i,j) = field_value[1];
      Ez(i,j) = field_value[2];
      Bx(i,j) = field_value[3];
      By(i,j) = field_value[4];
      Bz(i,j) = field_value[5];
    }
  }

  field->ComputeNormalizationFactor();

  // Output the data.
  x_field *= electronic_units.UNIT_LENGTH;
  y_field *= electronic_units.UNIT_LENGTH;

  x_field.save("x_field_qgauss.txt", arma::raw_ascii);
  y_field.save("y_field_qgauss.txt", arma::raw_ascii);
  Ex.save("QuasiGaussianField_Ex.txt", arma::raw_ascii);
  Ey.save("QuasiGaussianField_Ey.txt", arma::raw_ascii);
  Ez.save("QuasiGaussianField_Ez.txt", arma::raw_ascii);
  Bx.save("QuasiGaussianField_Bx.txt", arma::raw_ascii);
  By.save("QuasiGaussianField_By.txt", arma::raw_ascii);
  Bz.save("QuasiGaussianField_Bz.txt", arma::raw_ascii);

  // Compute the field in time at the focal spot.
  size_plot = 200;
  double tmax    = 5.0*lambda;
  auto   t_data  = arma::linspace<arma::colvec>(-tmax,tmax,size_plot);
  arma::colvec t_field(size_plot);

  for (uint i=0; i<size_plot; i++)
  {
    t_field(i) = field->ComputeFieldComponents(t_data(i), 0.0,0.0,0.0)[2];
  }

  t_data *= electronic_units.UNIT_TIME;
  t_data.save("t_data_qgauss.txt", arma::raw_ascii);
  t_field.save("t_field_qgauss.txt", arma::raw_ascii);
}

class DipoleQuasiGaussianQEDTest : public testing::Test
{
public:
  DipoleQuasiGaussianQEDTest()
  : UNIT_LENGTH(3.86159e-13)
  , UNIT_ENERGY(9.1093829140e-31*299792458.0*299792458.0)
  , UNIT_TIME(1.2880885e-21)
  , lambda(800e-9/UNIT_LENGTH)
  , omega(2.0*cst::pi<double>()/lambda)
  , pulse_duration(0.1*omega)
  , xmax(lambda)
  , energy(15.0/UNIT_ENERGY)
  , field(omega,pulse_duration,energy)
  {
    std::cout << "Energy [in QED Units]: " << energy << std::endl;
  }

protected:

  virtual void SetUp()
  {}

  const double UNIT_LENGTH;
  const double UNIT_ENERGY;
  const double UNIT_TIME;
  double lambda;
  double omega;

  double pulse_duration;
  double xmax;
  double energy;

  DipoleQuasiGaussian  field;
};

TEST_F(DipoleQuasiGaussianQEDTest, q_gauss)
{
  // Define the mesh of the plot.
  uint size_plot = 100;
  auto x_field = arma::linspace<arma::colvec>(-xmax,xmax, size_plot);
  auto y_field = arma::linspace<arma::colvec>(-xmax,xmax, size_plot);
  arma::mat field_values(size_plot,size_plot);

  // Compute the field values.
  double omega_0 = 2.0*constants::math::pi*constants::physics::c/800.0e-9;
  double mellotron_to_qed_factor = std::sqrt(4.0*constants::math::pi*constants::physics::alpha)*constants::physics::electron_mass*std::pow(constants::physics::c,2)/(omega_0*constants::physics::hbar);

  for (uint i=0; i<size_plot; i++)
  {
    for (uint j=0; j<size_plot; j++)
    {
      field_values(i,j) = mellotron_to_qed_factor*field.ComputeFieldComponents(0.0*lambda,x_field[i],y_field[j],0.0)[2];
    }
  }

  field.ComputeNormalizationFactor();

  // Output the data.
  x_field *= UNIT_LENGTH;
  y_field *= UNIT_LENGTH;

  x_field.save("x_field_qgauss_qed.txt", arma::raw_ascii);
  y_field.save("y_field_qgauss_qed.txt", arma::raw_ascii);
  field_values.save("QuasiGaussianField_qed.txt", arma::raw_ascii);

  // Compute the field in time at the focal spot.
  size_plot = 200;
  double tmax    = 5.0*lambda;
  auto   t_data  = arma::linspace<arma::colvec>(-tmax,tmax,size_plot);
  arma::colvec t_field(size_plot);

  for (uint i=0; i<size_plot; i++)
  {
    t_field(i) = mellotron_to_qed_factor*field.ComputeFieldComponents(t_data(i), 0.0,0.0,0.0)[2];
  }

  t_data *= UNIT_TIME;
  t_data.save("t_data_qgauss_qed.txt", arma::raw_ascii);
  t_field.save("t_field_qgauss_qed.txt", arma::raw_ascii);
}

GTEST_API_ int main(int argc, char **argv)
{
  H5open();
  printf("Running main() UnitTest-Particle.cpp.\n");
  testing::InitGoogleTest(&argc, argv);
  auto result =  RUN_ALL_TESTS();
  H5close();

  return result;
}


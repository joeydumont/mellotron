/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont                  <joey.dumont@gmail.com>               *
 * \since 2016-10-19                                                          *
 *                                                                            *
 * Unit test of fields/SalaminTightlyFocused in MELLOTRON. We simply evaluate *
 * and plot the field in a plane and compare with the original article.       *
 * --------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <armadillo>
#include <mellotron>

using namespace mellotron;

class SalaminTightlyFocusedTest : public testing::Test
{
public:
  SalaminTightlyFocusedTest()
  : lambda(800.0e-9)
  , omega_0(2.0*constants::math::pi*constants::physics::c/lambda)
  , electron_units(omega_0)
  , w0(0.272*lambda/electron_units.UNIT_LENGTH)
  , L(14.4*lambda/electron_units.UNIT_LENGTH)
  , xmax(1.5*lambda/electron_units.UNIT_LENGTH)
  , energy(13.0/electron_units.UNIT_ENERGY)
  , field(lambda/electron_units.UNIT_LENGTH,w0,L,energy)
  {
    lambda /= electron_units.UNIT_LENGTH;

    std::cout << "Norm factor: " << field.norm_factor << "  " <<  1.0/field.norm_factor << std::endl;
    std::cout << "Energy [in Mellotron Units]: " << energy << std::endl;
    std::cout << "Unit intensity: " << electron_units.UNIT_E_INTENSITY << std::endl;
  }

protected:

  virtual void SetUp()
  {}

  double lambda;
  const double omega_0;
  MellotronUnits                  electron_units;
  const double w0;
  const double L;
  const double xmax;
  const double energy;

  SalaminTightlyFocusedLinear     field;
};

TEST_F(SalaminTightlyFocusedTest, Linear)
{
  // Define the mesh of the plot.
  uint size_plot = 100;
  auto x_field = arma::linspace<arma::colvec>(-xmax,xmax, size_plot);
  auto y_field = arma::linspace<arma::colvec>(-xmax,xmax, 2*size_plot);
  arma::mat Ex(size_plot,2*size_plot);
  auto      Ey =Ex;
  auto      Ez =Ex;
  auto      Bx =Ex;
  auto      By =Ex;
  auto      Bz =Ex;

  // Compute the field values.
  for (uint i=0; i<size_plot; i++)
  {
    for (uint j=0; j<2*size_plot; j++)
    {
      auto field_vector = field.ComputeFieldComponents(lambda/4.0,x_field[i],y_field[j],lambda/2.0);
      Ex(i,j) = field_vector[0];
      Ey(i,j) = field_vector[1];
      Ez(i,j) = field_vector[2];
      Bx(i,j) = field_vector[3];
      By(i,j) = field_vector[4];
      Bz(i,j) = field_vector[5];
    }
  }

  // Compare to test data.
  arma::mat Ex_ref; Ex_ref.load("test_data/SalaminField_Ex.txt", arma::raw_ascii);
  arma::mat Ey_ref; Ey_ref.load("test_data/SalaminField_Ey.txt", arma::raw_ascii);
  arma::mat Ez_ref; Ez_ref.load("test_data/SalaminField_Ez.txt", arma::raw_ascii);
  arma::mat Bx_ref; Bx_ref.load("test_data/SalaminField_Bx.txt", arma::raw_ascii);
  arma::mat By_ref; By_ref.load("test_data/SalaminField_By.txt", arma::raw_ascii);
  arma::mat Bz_ref; Bz_ref.load("test_data/SalaminField_Bz.txt", arma::raw_ascii);

  for (uint i=0; i<size_plot; i++)
  {
    for (uint j=0; j<2*size_plot; j++)
    {
      EXPECT_NEAR(Ex(i,j), Ex_ref(i,j), 1.0e-3);
      EXPECT_NEAR(Ey(i,j), Ey_ref(i,j), 1.0e-3);
      EXPECT_NEAR(Ez(i,j), Ez_ref(i,j), 1.0e-3);
      EXPECT_NEAR(Bx(i,j), Bx_ref(i,j), 1.0e-3);
      EXPECT_NEAR(By(i,j), By_ref(i,j), 1.0e-3);
      EXPECT_NEAR(Bz(i,j), Bz_ref(i,j), 1.0e-3);
    }
  }

  // Output the data.
  x_field *= electron_units.UNIT_LENGTH;
  y_field *= electron_units.UNIT_LENGTH;

  x_field.save("x_field_salamin.txt", arma::raw_ascii);
  y_field.save("y_field_salamin.txt", arma::raw_ascii);
  Ex.save("SalaminField_Ex.txt", arma::raw_ascii);
  Ey.save("SalaminField_Ey.txt", arma::raw_ascii);
  Ez.save("SalaminField_Ez.txt", arma::raw_ascii);
  Bx.save("SalaminField_Bx.txt", arma::raw_ascii);
  By.save("SalaminField_By.txt", arma::raw_ascii);
  Bz.save("SalaminField_Bz.txt", arma::raw_ascii);

  // Compute intensity as a function of time.
  // I = 0.5c*epsilon_0*E^2.
  uint time_steps = 300;
  arma::vec intensityTimes  = arma::linspace<arma::vec>(-100e-15,100e-15,time_steps);
  arma::vec intensityValues(time_steps);

  for (uint i = 0; i < time_steps; i++)
  {
    auto field_vector  = field.ComputeFieldComponents(intensityTimes[i]*omega_0, 0.0,0.0,0.0);
    intensityValues[i] = 0.5*constants::physics::c*constants::physics::epsilon_0*std::pow(electron_units.UNIT_E_FIELD,2)*(field_vector[0]*field_vector[0]+field_vector[1]*field_vector[1]+field_vector[2]*field_vector[2])*1e-4;
  }

  intensityTimes.save("SalaminTimeIe.txt", arma::raw_ascii);
  intensityValues.save("SalaminTimeIe_time.txt", arma::raw_ascii);
}

class SalaminTightlyFocusedNormalizationTest : public testing::Test
{
public:
  SalaminTightlyFocusedNormalizationTest()
  : lambda(800e-9)
  , omega_0(2.0*constants::math::pi*constants::physics::c/lambda)
  , electron_units(omega_0)
  , w0(0.7*lambda/electron_units.UNIT_LENGTH)
  , L(0.8*lambda/electron_units.UNIT_LENGTH)
  , xmax(1.5*lambda/electron_units.UNIT_LENGTH)
  , energy(13.0/electron_units.UNIT_ENERGY)
  , field_norm_test(lambda/electron_units.UNIT_LENGTH,w0,L,1.0)
  , field(lambda/electron_units.UNIT_LENGTH,w0,L,field_norm_test.norm_factor,energy)
  {
    lambda /= electron_units.UNIT_LENGTH;

    std::cout << "Energy [in Mellotron Units]: " << energy << std::endl;
  }

protected:

  virtual void SetUp()
  {}

  double lambda;
  const double omega_0;
  MellotronUnits                  electron_units;
  const double w0;
  const double L;
  const double xmax;
  const double energy;

  SalaminTightlyFocusedLinear     field_norm_test;
  SalaminTightlyFocusedLinear     field;
};

TEST_F(SalaminTightlyFocusedNormalizationTest, Linear)
{
  // Define the mesh of the plot.
  uint size_plot = 200  ;
  auto x_field = arma::linspace<arma::colvec>(-xmax,xmax, size_plot);
  auto y_field = arma::linspace<arma::colvec>(-xmax,xmax, size_plot);
  arma::mat field_values(size_plot,size_plot);
  arma::mat field_values_norm(size_plot,size_plot);

  // Compute the field values.
  for (uint i=0; i<size_plot; i++)
  {
    for (uint j=0; j<size_plot; j++)
    {
      field_values(i,j)      = field.ComputeFieldComponents(0.0,x_field[i],y_field[j],lambda/(2.0*electron_units.UNIT_LENGTH))[0];
      field_values_norm(i,j) = field_norm_test.ComputeFieldComponents(0.0,x_field[i],y_field[j],lambda/(2.0*electron_units.UNIT_LENGTH))[0];

      //std::cout << field_values(i,j) << "\t" << std::sqrt(energy)*field_values_norm(i,j) << std::endl;

      EXPECT_NEAR(field_values(i,j),std::sqrt(energy)*field_values_norm(i,j),1.0e-5);
    }
  }

  // Output the data.
  x_field *= electron_units.UNIT_LENGTH;
  y_field *= electron_units.UNIT_LENGTH;

  //x_field.save("x_field_salamin.txt", arma::raw_ascii);
  //y_field.save("y_field_salamin.txt", arma::raw_ascii);
  field_values.save("SalaminField.txt", arma::raw_ascii);
}

class SalaminTightlyFocusedQEDTest : public testing::Test
{
public:
  SalaminTightlyFocusedQEDTest()
  : UNIT_LENGTH(3.86159e-13)
  , UNIT_ENERGY(9.1093829140e-31*299792458.0*299792458.0)
  , UNIT_TIME(1.2880885e-21)
  , lambda(800e-9/UNIT_LENGTH)
  , omega(2.0*cst::pi<double>()/lambda)
  , w0(0.7*lambda)
  , L(0.8*lambda)
  , xmax(1.5*lambda)
  , energy(15.0/UNIT_ENERGY)
  , field(lambda,w0,L,energy)
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
  const double omega;
  const double w0;
  const double L;
  const double xmax;
  const double energy;

  SalaminTightlyFocusedLinear     field;
};

TEST_F(SalaminTightlyFocusedQEDTest, Linear)
{
  // Define the mesh of the plot.
  uint size_plot = 100;
  auto x_field = arma::linspace<arma::colvec>(-xmax,xmax, size_plot);
  auto y_field = arma::linspace<arma::colvec>(-xmax,xmax, size_plot);
  arma::mat field_values(size_plot,size_plot);

  // Compute the field values.
  for (uint i=0; i<size_plot; i++)
  {
    for (uint j=0; j<size_plot; j++)
    {
      field_values(i,j) = field.ComputeFieldComponents(lambda/4.0,x_field[i],y_field[j],lambda/2.0)[0];
    }
  }

  field.ComputeNormalizationFactor();
  std::cout << "New norm factor should be unity: " << field.norm_factor << std::endl;

  // Output the data.
  x_field *= UNIT_LENGTH;
  y_field *= UNIT_LENGTH;

  //x_field.save("x_field_salamin.txt", arma::raw_ascii);
  //y_field.save("y_field_salamin.txt", arma::raw_ascii);
  field_values.save("SalaminField_qed.txt", arma::raw_ascii);
}

GTEST_API_ int main(int argc, char **argv)
{
  H5open();
  printf("Running main() UnitTest-SalaminTightlyFocused.cpp.\n");
  testing::InitGoogleTest(&argc, argv);
  auto result =  RUN_ALL_TESTS();
  H5close();

  return result;
}

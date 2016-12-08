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
  : lambda(800e-9)
  , omega_0(2.0*constants::math::pi*constants::physics::c/lambda)
  , electron_units(omega_0)
  , w0(0.7*lambda/electron_units.UNIT_LENGTH)
  , L(0.8*lambda/electron_units.UNIT_LENGTH)
  , xmax(1.5*lambda/electron_units.UNIT_LENGTH)
  , energy(15.0/electron_units.UNIT_ENERGY)
  , field(lambda/electron_units.UNIT_LENGTH,w0,L,energy)
  {
    lambda /= electron_units.UNIT_LENGTH;
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

  // Output the data.
  x_field *= electron_units.UNIT_LENGTH;
  y_field *= electron_units.UNIT_LENGTH;

  x_field.save("x_field_salamin.txt", arma::raw_ascii);
  y_field.save("y_field_salamin.txt", arma::raw_ascii);
  field_values.save("SalaminField.txt", arma::raw_ascii);
}

GTEST_API_ int main(int argc, char **argv)
{
  H5open();
  printf("Running main() UnitTest-Particle.cpp.\n");
  testing::InitGoogleTest(&argc, argv);
//  MPI_Init(&argc,&argv);
  auto result =  RUN_ALL_TESTS();
  //MPI_Finalize();
  H5close();

  return result;
}

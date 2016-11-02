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

namespace cst = boost::math::constants;

class DipoleQuasiGaussianTest : public testing::Test
{
public:
  DipoleQuasiGaussianTest()
  : UNIT_LENGTH(3.86159e-13)
  , UNIT_ENERGY(9.1093829140e-31*299792458.0*299792458.0)
  , UNIT_TIME(1.2880885e-21)
  , lambda(800e-9/UNIT_LENGTH)
  , omega(2.0*cst::pi<double>()/lambda)
  , pulse_duration(0.1*omega)
  , xmax(lambda)
  , energy(15.0/UNIT_ENERGY)
  , field(omega,pulse_duration,energy)
  {}

protected:

  virtual void SetUp()
  {}

  const double UNIT_LENGTH;
  const double UNIT_ENERGY;
  const double UNIT_TIME;
  const double lambda;
  const double omega;
  const double pulse_duration;
  const double xmax;
  const double energy;

  DipoleQuasiGaussian field;
};

TEST_F(DipoleQuasiGaussianTest, q_gauss)
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
      field_values(i,j) = field.ComputeFieldComponents(0.5*lambda,x_field[i],y_field[j],0.0)[2];
    }
  }

  field.ComputeNormalizationFactor();


  // Output the data.
  x_field *= UNIT_LENGTH;
  y_field *= UNIT_LENGTH;

  x_field.save("x_field_qgauss.txt", arma::raw_ascii);
  y_field.save("y_field_qgauss.txt", arma::raw_ascii);
  field_values.save("QuasiGaussianField.txt", arma::raw_ascii);

  // Compute the field in time at the focal spot.
  size_plot = 200;
  double tmax    = 5.0*lambda;
  auto   t_data  = arma::linspace<arma::colvec>(-tmax,tmax,size_plot);
  arma::colvec t_field(size_plot);

  for (uint i=0; i<size_plot; i++)
  {
    t_field(i) = field.ComputeFieldComponents(t_data(i), 0.0,0.0,0.0)[2];
  }

  t_data *= UNIT_TIME;
  t_data.save("t_data_qgauss.txt", arma::raw_ascii);
  t_field.save("t_field_qgauss.txt", arma::raw_ascii);
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

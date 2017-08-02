/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont                  <joey.dumont@gmail.com>               *
 * \since 2016-10-04                                                          *
 *                                                                            *
 * Unit test of Particle in MELLOTRON. The particle is instantiated as an     *
 * electron with original position at the origin and zero velocity.           *
 * --------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <armadillo>
#include <mellotron>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

using namespace mellotron;

struct ConstantField
{
  // Constructor to manually set the components of the static field.
  ConstantField(double my_Ex, double my_Ey, double my_Ez, double my_Bx, double my_By, double my_Bz)
  : Ex(my_Ex)
  , Ey(my_Ey)
  , Ez(my_Ez)
  , Bx(my_Bx)
  , By(my_By)
  , Bz(my_Bz)
  {}

  // Default constructor sets an electrostatic field.
  ConstantField()
  : Ex(0.0)
  , Ey(0.0)
  , Ez(0.1)
  , Bx(0.0)
  , By(0.0)
  , Bz(0.0)
  {}

  // Return the field components in an array.
  std::array<double,6> ComputeFieldComponents(double t, double x, double y, double z) const
  {
    std::array<double,6> constant = {Ex,Ey,Ez,Bx,By,Bz};
    return constant;
  }

  double Ex,Ey,Ez,Bx,By,Bz;
};

// We declare a test fixture to test a specific instance
// of Particle.
// This tests the standard Particle and ParticleObserver classes and compares
// their output to analytical solutions.
class ParticleTest : public testing::Test
{
public:
  ParticleTest()
  : charge(2.0)
  , mass(3.0)
  , electron_units(1.0)
  , electron(charge,mass,field,electron_units)
  , electron_obs(electron)
  {}

  const double charge;
  const double mass;

  ConstantField                                     field;
  MellotronUnits                                    electron_units;
  Particle<ConstantField>                           electron;
  ParticleObserver<ConstantField>                   electron_obs;

protected:

  virtual void SetUp()
  {
  }

};

// We test that we compute the proper coordinates.
// Analytical solution in an electrostatic field.
TEST_F(ParticleTest, TestIntegrationElectrostatic)
{
  // Define the initial conditions.
  arma::colvec::fixed<8> x           = arma::zeros<arma::colvec>(8);
  double x_init                      = 0.0;
  double y_init                      = 0.0;
  double z_init                      = 0.0;
  double px_init                     = 0.0;
  double py_init                     = 0.0;
  double pz_init                     = 0.0;

  // Times at which we output the data.
  unsigned int           size_time   = 2500;
  arma::colvec           times       = arma::linspace<arma::colvec>(0.0,10.0,size_time);

  // Set the initial conditions.
  electron.SetInitConditions(x,x_init,y_init,z_init,px_init,py_init,pz_init,times[0]);

  std::cout << "Initial conditions" << std::endl;
  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << std::endl;

  size_t steps = integrate_times(make_dense_output(1.0e-8,1.0e-8, runge_kutta_dopri5<arma::colvec::fixed<8> >() ),
                                 std::ref(electron),
                                 x,
                                 boost::begin(times),
                                 boost::end(times),
                                 0.01,
                                 std::ref(electron_obs));
  std::cout << steps << std::endl << std::endl;

  std::cout << "Final vector" << std::endl;
  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << std::endl;

  // Comparison between recorded data and analytical solution.
  for (uint i=0; i<size_time; i++)
  {
    EXPECT_NEAR(electron_obs.momentum(0,i), charge*field.Ex*electron_obs.times[i]+px_init*field.Ex, 1.0e-6);
    EXPECT_NEAR(electron_obs.momentum(1,i), charge*field.Ey*electron_obs.times[i]+py_init*field.Ey, 1.0e-6);
    EXPECT_NEAR(electron_obs.momentum(2,i), charge*field.Ez*electron_obs.times[i]+pz_init*field.Ez, 1.0e-6);

    double field_norm_sq = std::pow(arma::norm(electron_obs.electric_field.col(i),2),2);
    EXPECT_NEAR(electron_obs.position(0,i), field.Ex/field_norm_sq*mass/charge*(electron_obs.gamma[i]-electron_obs.gamma[0])+x_init, 1.0e-6);
    EXPECT_NEAR(electron_obs.position(1,i), field.Ey/field_norm_sq*mass/charge*(electron_obs.gamma[i]-electron_obs.gamma[0])+y_init, 1.0e-6);
    EXPECT_NEAR(electron_obs.position(2,i), field.Ez/field_norm_sq*mass/charge*(electron_obs.gamma[i]-electron_obs.gamma[0])+z_init, 1.0e-6);
  }

  electron_obs.OutputData();
}

// We test that we compute the right solution for a uniform magnetostatic field.
TEST_F(ParticleTest, TestIntegrationMagnetostatic)
{
 // Define the initial conditions.
  arma::colvec::fixed<8> x           = arma::zeros<arma::colvec>(8);
  double x_init                      = 0.0;
  double y_init                      = 0.0;
  double z_init                      = 0.0;
  double px_init                     = 0.1;
  double py_init                     = 0.1;
  double pz_init                     = 0.0;

  // Times at which we output the data.
  unsigned int           size_time   = 250;
  arma::colvec           times       = arma::linspace<arma::colvec>(0.0,10.0,size_time);

  // Set the initial conditions.
  electron.SetInitConditions(x,x_init,y_init,z_init,px_init,py_init,pz_init,times[0]);

  // Set the field.
  field.Ez=0.0;
  field.Bz=1.0;
  double omega   = charge*field.Bz/x[4];
  double theta_0 = std::atan2(px_init,py_init);
  double p_perp  = std::sqrt(px_init*px_init+py_init*py_init);

  std::cout << "Initial conditions" << std::endl;
  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << std::endl;

  size_t steps = integrate_times(make_dense_output(1.0e-8,1.0e-8, runge_kutta_dopri5<arma::colvec::fixed<8> >() ),
                                 std::ref(electron),
                                 x,
                                 boost::begin(times),
                                 boost::end(times),
                                 0.01,
                                 std::ref(electron_obs));
  std::cout << steps << std::endl << std::endl;

  std::cout << "Final vector" << std::endl;
  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << std::endl;

  // Comparison between recorded data and analytical solution.
  for (uint i=0; i<size_time; i++)
  {
    EXPECT_NEAR(electron_obs.momentum(0,i), p_perp*std::sin(omega*electron_obs.times[i]+theta_0), 1.0e-5);
    EXPECT_NEAR(electron_obs.momentum(1,i), p_perp*std::cos(omega*electron_obs.times[i]+theta_0), 1.0e-5);
    EXPECT_NEAR(electron_obs.momentum(2,i), pz_init, 1.0e-6);

    EXPECT_NEAR(electron_obs.position(0,i), -p_perp/(omega*electron_obs.gamma[i]*mass)*(cos(omega*electron_obs.times[i]+theta_0)-cos(theta_0)) + x_init, 1.0e-5);
    EXPECT_NEAR(electron_obs.position(1,i),  p_perp/(omega*electron_obs.gamma[i]*mass)*(sin(omega*electron_obs.times[i]+theta_0)-sin(theta_0)) + y_init, 1.0e-5);
    EXPECT_NEAR(electron_obs.position(2,i),  pz_init/(electron_obs.gamma[i]*mass)*electron_obs.times[i]+z_init, 1.0e-5);
  }
}

// We declare a test fixture to test a specific instance
// of ParticleIonized. This tests ParticleIonized and the relevant
// ParticleObserverIonized. This test should generate a single
// HDF5 file, as the second particle (electron_below), is always below
// threshold and should therefore not output data.
class ParticleIonizedTest : public testing::Test
{
public:
  ParticleIonizedTest()
  : charge(2.0)
  , mass(3.0)
  , electron_units(1.0)
  , electron(charge,mass,field,electron_units,0.0)
  , electron_obs(electron)
  , electron_below(charge,mass,field,electron_units,10000.0)
  , electron_obs_below(electron_below)
  {}

  const double charge;
  const double mass;

  ConstantField                                     field;
  MellotronUnits                                    electron_units;
  ParticleIonized<ConstantField>                    electron;
  ParticleObserverIonized<ConstantField>            electron_obs;
  ParticleIonized<ConstantField>                    electron_below;
  ParticleObserverIonized<ConstantField>            electron_obs_below;

protected:

  virtual void SetUp()
  {
  }

};

// We test that we compute the proper coordinates.
// Should output 8503784799965847432.hdf5.
TEST_F(ParticleIonizedTest, TestIntegrationElectrostatic)
{
  // Define the initial conditions.
  arma::colvec::fixed<8> x           = arma::zeros<arma::colvec>(8);
  double x_init                      = 0.0;
  double y_init                      = 0.0;
  double z_init                      = 0.0;
  double px_init                     = 0.0;
  double py_init                     = 0.0;
  double pz_init                     = 0.0;

  // Times at which we output the data.
  unsigned int           size_time   = 2500;
  arma::colvec           times       = arma::linspace<arma::colvec>(0.0,10.0,size_time);

  // Set the initial conditions.
  electron.SetInitConditions(x,x_init,y_init,z_init,px_init,py_init,pz_init,times[0]);

  std::cout << "Initial conditions" << std::endl;
  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << std::endl;

  size_t steps = integrate_times(make_dense_output(1.0e-7,1.0e-7, runge_kutta_dopri5<arma::colvec::fixed<8> >() ),
                                 std::ref(electron),
                                 x,
                                 boost::begin(times),
                                 boost::end(times),
                                 0.01,
                                 std::ref(electron_obs));
  std::cout << steps << std::endl << std::endl;

  std::cout << "Final vector" << std::endl;
  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << std::endl;

  // Comparison between recorded data and analytical solution.
  for (uint i=0; i<size_time; i++)
  {
    EXPECT_NEAR(electron_obs.momentum(0,i), charge*field.Ex*electron_obs.times[i]+px_init*field.Ex, 1.0e-6);
    EXPECT_NEAR(electron_obs.momentum(1,i), charge*field.Ey*electron_obs.times[i]+py_init*field.Ey, 1.0e-6);
    EXPECT_NEAR(electron_obs.momentum(2,i), charge*field.Ez*electron_obs.times[i]+pz_init*field.Ez, 1.0e-6);

    double field_norm_sq = std::pow(arma::norm(electron_obs.electric_field.col(i),2),2);
    EXPECT_NEAR(electron_obs.position(0,i), field.Ex/field_norm_sq*mass/charge*(electron_obs.gamma[i]-electron_obs.gamma[0])+x_init, 1.0e-6);
    EXPECT_NEAR(electron_obs.position(1,i), field.Ey/field_norm_sq*mass/charge*(electron_obs.gamma[i]-electron_obs.gamma[0])+y_init, 1.0e-6);
    EXPECT_NEAR(electron_obs.position(2,i), field.Ez/field_norm_sq*mass/charge*(electron_obs.gamma[i]-electron_obs.gamma[0])+z_init, 1.0e-6);

  }

  electron_obs.OutputData();
}

// Test a below threshold integration with ParticleIonized. Should
// not output anything. The different initial condition makes sure
// that any output would be unique.
TEST_F(ParticleIonizedTest, TestIntegrationElectrostaticBelowThreshold)
{
  // Define the initial conditions.
  arma::colvec::fixed<8> x           = arma::zeros<arma::colvec>(8);
  double x_init                      = 0.1;
  double y_init                      = 0.0;
  double z_init                      = 0.0;
  double px_init                     = 0.0;
  double py_init                     = 0.0;
  double pz_init                     = 0.0;

  // Times at which we output the data.
  unsigned int           size_time   = 2500;
  arma::colvec           times       = arma::linspace<arma::colvec>(0.0,10.0,size_time);

  // Set the initial conditions.
  electron_below.SetInitConditions(x,x_init,y_init,z_init,px_init,py_init,pz_init,times[0]);

  std::cout << "Initial conditions" << std::endl;
  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << std::endl;

  size_t steps = integrate_times(make_dense_output(1.0e-7,1.0e-7, runge_kutta_dopri5<arma::colvec::fixed<8> >() ),
                                 std::ref(electron_below),
                                 x,
                                 boost::begin(times),
                                 boost::end(times),
                                 0.01,
                                 std::ref(electron_obs_below));
  std::cout << steps << std::endl << std::endl;

  std::cout << "Final vector" << std::endl;
  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << std::endl;

  // Comparison between recorded data and analytical solution.
  for (uint i=0; i<size_time; i++)
  {
    EXPECT_NEAR(electron_obs_below.momentum(0,i), px_init, 1.0e-6);
    EXPECT_NEAR(electron_obs_below.momentum(1,i), py_init, 1.0e-6);
    EXPECT_NEAR(electron_obs_below.momentum(2,i), pz_init, 1.0e-6);

    double field_norm_sq = std::pow(arma::norm(electron_obs_below.electric_field.col(i),2),2);
    EXPECT_NEAR(electron_obs_below.position(0,i), px_init*times[i]+x_init, 1.0e-6);
    EXPECT_NEAR(electron_obs_below.position(1,i), py_init*times[i]+y_init, 1.0e-6);
    EXPECT_NEAR(electron_obs_below.position(2,i), pz_init*times[i]+z_init, 1.0e-6);

  }

  electron_obs_below.OutputData();
}

// We test that we compute the right solution for a uniform magnetostatic field
// with the ParticleIonized class.
TEST_F(ParticleIonizedTest, TestIntegratinoMagnetostatic)
{
 // Define the initial conditions.
  arma::colvec::fixed<8> x           = arma::zeros<arma::colvec>(8);
  double x_init                      = 0.0;
  double y_init                      = 0.0;
  double z_init                      = 0.0;
  double px_init                     = 0.1;
  double py_init                     = 0.1;
  double pz_init                     = 0.0;

  // Times at which we output the data.
  unsigned int           size_time   = 250;
  arma::colvec           times       = arma::linspace<arma::colvec>(0.0,10.0,size_time);

  // Set the initial conditions.
  electron.SetInitConditions(x,x_init,y_init,z_init,px_init,py_init,pz_init,times[0]);

  // Set the field.
  field.Ez=0.0;
  field.Bz=1.0;
  double omega   = charge*field.Bz/x[4];
  double theta_0 = std::atan2(px_init,py_init);
  double p_perp  = std::sqrt(px_init*px_init+py_init*py_init);

  std::cout << "Initial conditions" << std::endl;
  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << std::endl;

  size_t steps = integrate_times(make_dense_output(1.0e-7,1.0e-7, runge_kutta_dopri5<arma::colvec::fixed<8> >() ),
                                 std::ref(electron),
                                 x,
                                 boost::begin(times),
                                 boost::end(times),
                                 0.01,
                                 std::ref(electron_obs));
  std::cout << steps << std::endl << std::endl;

  std::cout << "Final vector" << std::endl;
  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
  std::cout << std::endl;

  // Comparison between recorded data and analytical solution.
  for (uint i=0; i<size_time; i++)
  {
    EXPECT_NEAR(electron_obs.momentum(0,i), p_perp*std::sin(omega*electron_obs.times[i]+theta_0), 1.0e-5);
    EXPECT_NEAR(electron_obs.momentum(1,i), p_perp*std::cos(omega*electron_obs.times[i]+theta_0), 1.0e-5);
    EXPECT_NEAR(electron_obs.momentum(2,i), pz_init, 1.0e-6);

    EXPECT_NEAR(electron_obs.position(0,i), -p_perp/(omega*electron_obs.gamma[i]*mass)*(cos(omega*electron_obs.times[i]+theta_0)-cos(theta_0)) + x_init, 1.0e-5);
    EXPECT_NEAR(electron_obs.position(1,i),  p_perp/(omega*electron_obs.gamma[i]*mass)*(sin(omega*electron_obs.times[i]+theta_0)-sin(theta_0)) + y_init, 1.0e-5);
    EXPECT_NEAR(electron_obs.position(2,i),  pz_init/(electron_obs.gamma[i]*mass)*electron_obs.times[i]+z_init, 1.0e-5);
  }
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

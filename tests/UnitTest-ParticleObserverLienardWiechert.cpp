/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont                  <joey.dumont@gmail.com>               *
 * \since 2017-07-23                                                          *
 *                                                                            *
 * Unit test of ParticleObserverLineardWiechert in MELLOTRON. The particle    *
 * electron in an electric field. We test whether the class                   *
 * ParticleObserverLienardWiechert properly outputs the field. Later tests    *
 * will compare against a reference simulation.                               *
 * --------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <cmath>

#include <armadillo>
#include <mellotron>
#include <boost/numeric/odeint.hpp>


using namespace boost::numeric::odeint;

using namespace mellotron;

struct PlaneWave
{
  PlaneWave(double my_omega)
  : Ex(0.0)
  , Ey(0.0)
  , Ez(0.1)
  , Bx(0.0)
  , By(0.0)
  , Bz(0.0)
  , omega(my_omega)
  {}

  std::array<double,6> ComputeFieldComponents(double t, double x, double y, double z) const
  {
    std::array<double,6> constant = {Ex,Ey,Ez*std::cos(omega*t),Bx,By,Bz};
    return constant;
  }

  double Ex,Ey,Ez,Bx,By,Bz,omega;
};

// We declare a test fixture to test a specific instance
// of Particle.
class ParticleObserverTest : public testing::Test
{
public:
  ParticleObserverTest()
  : charge(2.0)
  , mass(1.0)
  , radius(1e4)
  , number_points_theta(50)
  , number_points_phi(50)
  , field(1.0)
  , electron_units(1.0)
  , electron(charge,mass,field,electron_units)
  , electron_obs(electron,radius,number_points_theta,number_points_phi,100)
  {}

  const double charge;
  const double mass;
  const double radius;
  const double number_points_theta;
  const double number_points_phi;

  PlaneWave                                         field;
  MellotronUnits                                    electron_units;
  Particle<PlaneWave>                               electron;
  ParticleObserverLienardWiechert<PlaneWave>        electron_obs;

protected:

  virtual void SetUp()
  {
  }

};

// We test that we compute the proper coordinates.
TEST_F(ParticleObserverTest, TestConstantAcceleration)
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
  unsigned int           size_time   = 500;
  arma::colvec           times       = arma::linspace<arma::colvec>(0.0,100.0,size_time);

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

  // Check whether some data was computed.
  std::ofstream field_file;
  std::ofstream time_file;
  field_file.open("field.txt");
  time_file.open("time.txt");
  for (unsigned int i=0; i<size_time; i++)
  {
    field_file << electron_obs.electric_field_lw[i][number_points_theta/2][0][0] << "\n";
    time_file  <<  times(i) << "\n";
  }
  field_file.close();
  time_file.close();

  electron_obs.OutputData();
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

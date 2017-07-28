/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont                  <joey.dumont@gmail.com>               *
 * \since 2017-.7-23                                                          *
 *                                                                            *
 * Unit test of ParticleObserverLineardWiechert in MELLOTRON. The particle    *
 * is an electron in a magnetic field. We test whether the class              *
 * ParticleObserverLienardWiechert properly computes and outputs the field.   *
 * --------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <algorithm>

#include <armadillo>
#include <mellotron>
#include <boost/numeric/odeint.hpp>


using namespace boost::numeric::odeint;

using namespace mellotron;

struct ConstantField
{
  ConstantField(double my_Ex, double my_Ey, double my_Ez, double my_Bx, double my_By, double my_Bz)
  : Ex(my_Ex)
  , Ey(my_Ey)
  , Ez(my_Ez)
  , Bx(my_Bx)
  , By(my_By)
  , Bz(my_Bz)
  {}

  ConstantField()
  : Ex(0.0)
  , Ey(0.0)
  , Ez(0.1)
  , Bx(0.0)
  , By(0.0)
  , Bz(0.0)
  {}

  std::array<double,6> ComputeFieldComponents(double t, double x, double y, double z) const
  {
    std::array<double,6> constant = {Ex,Ey,Ez,Bx,By,Bz};
    return constant;
  }

  double Ex,Ey,Ez,Bx,By,Bz;
};

// We declare a test fixture to test a specific instance
// of Particle.
class ParticleTest : public testing::Test
{
public:
  ParticleTest()
  : charge(2.0)
  , mass(3.0)
  , radius(1e4)
  , number_points_theta(50)
  , number_points_phi(50)
  , electron_units(1.0)
  , electron(charge,mass,field,electron_units)
  , electron_obs(electron,radius,number_points_theta,number_points_phi,100)
  {}

  const double charge;
  const double mass;
  const double radius;
  const uint   number_points_theta;
  const uint   number_points_phi;

  ConstantField                                     field;
  MellotronUnits                                    electron_units;
  Particle<ConstantField>                           electron;
  ParticleObserverLienardWiechert<ConstantField>    electron_obs;

protected:

  virtual void SetUp()
  {
  }

};


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
    arma::colvec position(3);
    position(0) = -p_perp/(omega*electron_obs.gamma[i]*mass)*(cos(omega*electron_obs.times[i]+theta_0)-cos(theta_0)) + x_init;
    position(1) =  p_perp/(omega*electron_obs.gamma[i]*mass)*(sin(omega*electron_obs.times[i]+theta_0)-sin(theta_0)) + y_init;
    position(2) =  pz_init/(electron_obs.gamma[i]*mass)*electron_obs.times[i]+z_init;

    arma::colvec momentum(3);
    momentum(0) = p_perp*std::sin(omega*electron_obs.times[i]+theta_0);
    momentum(1) = p_perp*std::cos(omega*electron_obs.times[i]+theta_0);
    momentum(2) = pz_init;

    arma::colvec momentump(3);
    momentump(0) =  p_perp*omega*std::cos(omega*electron_obs.times[i]+theta_0);
    momentump(1) = -p_perp*omega*std::sin(omega*electron_obs.times[i]+theta_0);
    momentump(2) = 0.0;

    EXPECT_NEAR(electron_obs.momentum(0,i), momentum(0), 1.0e-5);
    EXPECT_NEAR(electron_obs.momentum(1,i), momentum(1), 1.0e-5);
    EXPECT_NEAR(electron_obs.momentum(2,i), momentum(2), 1.0e-6);

    EXPECT_NEAR(electron_obs.position(0,i),  position(0), 1.0e-5);
    EXPECT_NEAR(electron_obs.position(1,i),  position(1), 1.0e-5);
    EXPECT_NEAR(electron_obs.position(2,i),  position(2), 1.0e-5);

    // Find the maximum of the LW electric field.
    double max_e_field = *(std::max_element(electron_obs.electric_field_lw.origin(),
                                          electron_obs.electric_field_lw.origin() + electron_obs.electric_field_lw.num_elements()));

    // Now compute the Liénard-Wiechert fields.
    for (uint j=0; j<electron_obs.number_points_theta; j++)
    {
      double theta = electron_obs.theta[j];
      for (uint k=0; k<electron_obs.number_points_phi; k++)
      {
        double phi = electron_obs.phi[k];

        // Compute the normal as simply r.
        arma::colvec normal(3);
        normal(0) = std::sin(theta)*std::cos(phi);
        normal(1) = std::sin(theta)*std::sin(phi);
        normal(2) = std::cos(theta);

        // We compute the prefactor.
        double prefactor = constants::physics::hbar/(constants::physics::electron_mass*std::pow(constants::physics::c,2))
                                     * constants::physics::alpha * charge / electron_obs.radius;
        double denominator = std::pow(1.0-arma::dot(normal,momentum)/(electron_obs.gamma[i]*mass),-3.0);

        // We compute the field.
        arma::colvec firstTerm = normal-momentum/(electron_obs.gamma[i]*mass);
        arma::colvec secondTerm= momentum/(electron_obs.gamma[i]*mass) - arma::dot(momentum,momentump)*momentum/std::pow(electron_obs.gamma[i]*mass,3);

        arma::colvec e_field_lw = prefactor*denominator*arma::cross(normal, arma::cross(firstTerm,secondTerm));

        EXPECT_NEAR(electron_obs.electric_field_lw[i][j][k][0]/max_e_field, e_field_lw(0)/max_e_field, 1.0e-4);
        EXPECT_NEAR(electron_obs.electric_field_lw[i][j][k][1]/max_e_field, e_field_lw(1)/max_e_field, 1.0e-4);
        EXPECT_NEAR(electron_obs.electric_field_lw[i][j][k][2]/max_e_field, e_field_lw(2)/max_e_field, 1.0e-4);
      }
    }
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

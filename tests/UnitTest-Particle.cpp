/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont                  <joey.dumont@gmail.com>               *
 * \since 2016-10-04                                                          *
 *                                                                            *
 * Unit test of Particle in MELLOTRON. The particle is instantiated as an     *
 * electron with original position at the origin and zero velocity.           *
 * --------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <armadillo>
#include <mellotron>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;

struct ConstantField
{
  std::array<double,6> ComputeFieldComponents(double t, double x, double y, double z)
  {
    std::array<double,6> constant = {0.0,0.0,1.0e-1,0.0,0.0,1.0e-1};
    return constant;
  }
};

// We declare a test fixture to test a specific instance
// of Particle.
class ParticleTest : public testing::Test
{
public:
  ParticleTest()
  : electron(1.0,1.0,field)
  , electron_obs(electron)
  {}

protected:

  virtual void SetUp()
  {
    //electron     = new Particle<ConstantField>(1.0,1.0,field);
    //electron_obs = ParticleObserver<ConstantField>(*electron);
    //electron = new Particle<ConstantField>(1.0,1.0,field,std::string("ll_first_term"));
  }

  ConstantField                       field;
  Particle<ConstantField>             electron;
  ParticleObserver<ConstantField>     electron_obs;

};

// We test that we compute the proper coordinates.
TEST_F(ParticleTest, TestIntegration)
{
  // Define the initial conditions.
  arma::colvec::fixed<8> x           = arma::zeros<arma::colvec>(8);
  arma::colvec           times       = arma::linspace<arma::colvec>(0.0,10.0,50);
  x[4] = 1.0;
  size_t steps = integrate_times(make_dense_output(1.0e-6,1.0e-6, runge_kutta_dopri5<arma::colvec::fixed<8> >() ),
                                 std::ref(electron),
                                 x,
                                 boost::begin(times),
                                 boost::end(times),
                                 0.1,
                                 std::ref(electron_obs));
  std::cout << steps << std::endl;

  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }

  electron_obs.OutputData();
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
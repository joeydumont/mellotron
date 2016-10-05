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
  ParticleTest() {}

protected:

  virtual void SetUp()
  {
    electron = new Particle<ConstantField>(1.0,1.0,field);
    //electron = new Particle<ConstantField>(1.0,1.0,field,std::string("ll_first_term"));
  }

  ConstantField               field;
  Particle<ConstantField>  *  electron;

};

// We test that we compute the proper coordinates.
TEST_F(ParticleTest, TestIntegration)
{
  // Define the initial conditions.
  arma::colvec::fixed<8> x = arma::zeros<arma::colvec>(8);
  x[4] = 1.0;
  size_t steps = integrate(*electron, x, 0.0, 10.0, 0.1);
  std::cout << steps << std::endl;

  for (uint i=0; i<8; i++)
  {
    std::cout << x[i] << std::endl;
  }
}

GTEST_API_ int main(int argc, char **argv)
{
  printf("Running main() UnitTest-Particle.cpp.\n");
  testing::InitGoogleTest(&argc, argv);
//  MPI_Init(&argc,&argv);
  auto result =  RUN_ALL_TESTS();
  //MPI_Finalize();

  return result;
}
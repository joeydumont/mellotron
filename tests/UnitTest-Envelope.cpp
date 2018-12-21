/*! ------------------------------------------------------------------------- *
 * \author Francois Fillion-Gourdeau  <francois.fillion@emt.inrs.ca>          *
 * \since 2018-12-21                                                          *
 *                                                                            *
 * Unit test of Envelope in MELLOTRON.                                        *
 * --------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <mellotron>



using namespace mellotron;

TEST(Envelope, NoEnvelope)
{
  NoEnvelope env;
  EXPECT_NEAR(1.0, env.Value(0.0),1.0e-10);
}
TEST(Envelope, EnvelopeHann)
{
  EnvelopeHann env(-0.5,1.0);
  EXPECT_NEAR(1.0, env.Value(0.0),1.0e-10);
  EXPECT_NEAR(0.0, env.Value(0.5),1.0e-10);
  EXPECT_NEAR(0.0, env.Value(-0.5),1.0e-10);
  EXPECT_NEAR(0.5, env.Value(0.25),1.0e-10);
  EXPECT_NEAR(0.5, env.Value(-0.25),1.0e-10);
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

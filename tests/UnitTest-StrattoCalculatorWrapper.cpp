/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont                  <joey.dumont@gmail.com>               *
 * \since 2016-10-19                                                          *
 *                                                                            *
 * Unit test of simulations/StrattoCalculatorWrapper in MELLOTRON. We simply  *
 * evaluate and plot the field in the focal plane and compare to results      *
 * obtained directly with the StrattoCalculator.                              *
 * --------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "../simulations/StrattoCalculatorWrapper.hpp"
#include <armadillo>
#include <mellotron>
#include <meshpi>
#include <strattocalculator>

using namespace mellotron;
using namespace MeshPI;
using namespace StrattoCalculator;

class StrattoCalculatorWrapperTest : public testing::Test
{
public:
  StrattoCalculatorWrapperTest()
  : lambda(800.0e-9)
  , omega_0(2.0*constants::math::pi*constants::physics::c/lambda)
  , electron_units(omega_0)
  , beam_width(0.075/electron_units.UNIT_LENGTH)
  , focal_length(0.04375/electron_units.UNIT_LENGTH)
  , xmax(2.5*lambda/electron_units.UNIT_LENGTH)
  {

    // Create the mesh of the parabola.
    std::array<double,2> min_parabola                      = {0.0,0.0};
    std::array<double,2> max_parabola                      = {0.0875/electron_units.UNIT_LENGTH,6.28318530718};
    std::array<unsigned int, 2> intervals_per_dim_parabola = {25,25};
    coord_sys       = new CoordinateSystemModel<2,Cylindrical>(x12,space,qed_units);
    domain_parabola                       = new DomainSerial<2>(
                                                                min_parabola,
                                                                max_parabola,
                                                                intervals_per_dim_parabola,
                                                                coord_sys
                                                                );

    std::array<unsigned int, 2> num_points_parabola        = {5,5};
    std::array<bool,2> boundary                            = {true,false};
    mesh_parabola                                          = new SurfaceMesh<2>(
        domain_parabola,
        num_points_parabola,
        MeshPI::gausslegendre_mesh,
        boundary
    );

    surface                    = new EmittingSurfaceParabola<2>(
        mesh_parabola,
        focal_length
        );
        surface->ComputeSurface();
        surface->ComputeNormal();

    // Change spectrum parameters to appropriate units.
    double energy       = 13.0/electron_units.UNIT_ENERGY;
    double lambda_c     = lambda/electron_units.UNIT_LENGTH;
    double lambda_min   = 720.0e-9/electron_units.UNIT_LENGTH;
    double lambda_max   = 880.0e-9/electron_units.UNIT_LENGTH;
    double delta_lambda = 60.0e-9/electron_units.UNIT_LENGTH;

    // Create the spectrum.
    int number_components = 25;
    std::vector<double> beam_frequencies, beam_freq_widths, beam_freq_ratios;
    beam_frequencies.push_back(lambda_c);
    beam_freq_widths.push_back(delta_lambda);
    beam_freq_ratios.push_back(1.0);
    std::vector<double> phase(number_components,0.0);
    spectrum_incident = new SpectrumGaussianWavelength(
        number_components,
        lambda_min,
        lambda_max,
        phase,
        beam_frequencies,
        beam_freq_widths,
        beam_freq_ratios,
        energy,
        8
    );
    for (int i=0; i<number_components;i++)
    {
        spectrum_incident->SetPhase(i, -2.0*spectrum_incident->GetOmega(i)*focal_length);
    }

    // Create the beam model.
    beam = new TEM00Mode(spectrum_incident,beam_width);

    // Evaluate the beam on the mirror.
    incident_field = new SurfaceEMFieldManyOnTheFly<SurfaceEMFieldGeneral,2,1>(
        mesh_parabola,
        spectrum_incident,
        surface,
        beam
    );

    // Declare the integrator and convert it to the form the MELLOTRON expects.
    integrator          = new StrattonChuIntegrator::MeshlessAxialSurfGenField(surface,incident_field);
    meshless_time_field = new TemporalEMFieldMeshless<StrattonChuIntegrator::MeshlessAxialSurfGenField>(integrator);
    field               = new StrattoCalculatorWrapper<StrattonChuIntegrator::MeshlessAxialSurfGenField>(*meshless_time_field);
  }

  ~StrattoCalculatorWrapperTest()
  {
    delete  coord_sys;
    delete  domain_parabola;
    delete  mesh_parabola;
    delete  surface;
    //delete  spectrum_incident;
    delete  beam;
    delete  incident_field;
    delete  integrator;
    delete  meshless_time_field;
    delete  field;
  }

protected:

  virtual void SetUp()
  {}

        double          lambda;
  const double          omega_0;
        MellotronUnits  electron_units;
  const double          beam_width,focal_length;
  const double          xmax;
        double          lambda_c;
        double          energy;

  CoordinateSystemModel<2,Cylindrical>                                        *  coord_sys;
  DomainSerial<2>                                                             *  domain_parabola;
  SurfaceMesh<2>                                                              *  mesh_parabola;
  EmittingSurfaceAxialSym<2>                                                  *  surface;
  SpectrumGaussianWavelength                                                  *  spectrum_incident;
  TEM00Mode                                                                   *  beam;
  SurfaceEMFieldManyOnTheFly<SurfaceEMFieldGeneral,2,1>                       *  incident_field;
  StrattonChuIntegrator::MeshlessAxialSurfGenField                            *  integrator;
  TemporalEMFieldMeshless<StrattonChuIntegrator::MeshlessAxialSurfGenField>   *  meshless_time_field;
  StrattoCalculatorWrapper<StrattonChuIntegrator::MeshlessAxialSurfGenField>  *  field;
};

TEST_F(StrattoCalculatorWrapperTest, Linear)
{
  // Define the mesh of the plot.
  uint size_plot = 50;
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
      auto field_vector = field->ComputeFieldComponents(0.0,x_field[i],y_field[j],0.0);
      Ex(i,j) = field_vector[0];
      Ey(i,j) = field_vector[1];
      Ez(i,j) = field_vector[2];
      Bx(i,j) = field_vector[3];
      By(i,j) = field_vector[4];
      Bz(i,j) = field_vector[5];
    }
  }

  // Compute Ex at z=0 in as a function of time.
  uint size_time = 1000;
  double tmax    = -200.0;
  auto t_field   = arma::linspace(-tmax,tmax,size_time);
  arma::mat Ex_t(size_time,4);

  for (uint i=0; i<size_time; i++)
  {
    Ex_t(i,0) = t_field(i);
    auto field_time = field->ComputeFieldComponents(t_field(i),0.0,0.0,0.0);
    Ex_t(i,arma::span(1,3)) = arma::rowvec({field_time[0],field_time[1],field_time[2]});
  }

  // Output the data.
  x_field *= electron_units.UNIT_LENGTH;
  y_field *= electron_units.UNIT_LENGTH;

  x_field.save("x_field_stratto.txt", arma::raw_ascii);
  y_field.save("y_field_stratto.txt", arma::raw_ascii);
  Ex.save("StrattoField_Ex.txt", arma::raw_ascii);
  Ey.save("StrattoField_Ey.txt", arma::raw_ascii);
  Ez.save("StrattoField_Ez.txt", arma::raw_ascii);
  Bx.save("StrattoField_Bx.txt", arma::raw_ascii);
  By.save("StrattoField_By.txt", arma::raw_ascii);
  Bz.save("StrattoField_Bz.txt", arma::raw_ascii);

  Ex_t.save("StrattoField_Ex_time.txt", arma::raw_ascii);

  // Compare to test data.
  arma::mat Ex_ref; Ex_ref.load("test_data/StrattoField_Ex.txt", arma::raw_ascii);
  arma::mat Ey_ref; Ey_ref.load("test_data/StrattoField_Ey.txt", arma::raw_ascii);
  arma::mat Ez_ref; Ez_ref.load("test_data/StrattoField_Ez.txt", arma::raw_ascii);
  arma::mat Bx_ref; Bx_ref.load("test_data/StrattoField_Bx.txt", arma::raw_ascii);
  arma::mat By_ref; By_ref.load("test_data/StrattoField_By.txt", arma::raw_ascii);
  arma::mat Bz_ref; Bz_ref.load("test_data/StrattoField_Bz.txt", arma::raw_ascii);

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
}

GTEST_API_ int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  H5open();
  printf("Running main() UnitTest-StrattoCalculatorWrapper.cpp.\n");
  testing::InitGoogleTest(&argc, argv);
  auto result =  RUN_ALL_TESTS();
  H5close();
  MPI_Finalize();

  return result;
}

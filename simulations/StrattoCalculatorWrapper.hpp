/*! ------------------------------------------------------------------------- *
 * \author Joey Dumont  										              *
 * \since 2017-05-31                                                          *
 *                                                                            *
 * Header-only wrapper to use the field computation abilities of the          *
 * StrattoCalculator from the MELLOTRON.                                      *
 * --------------------------------------------------------------------------*/

#include <armadillo>
#include <meshpi>
#include <strattocalculator>
#include <mellotron>

using namespace MeshPI;
using namespace StrattoCalculator;

template <class FieldRepresentation>
class StrattoCalculatorWrapper
{
public:
    /// Sets the field representation.
    StrattoCalculatorWrapper(TemporalEMFieldMeshless<FieldRepresentation> & my_field_rep)
    : field_rep(my_field_rep)
    {}

    /// Computes the field in Cartesian coordinates, as the MELLOTRON expects.
    std::array<double,6> ComputeFieldComponents(double t, double x, double y, double z)
    {
        // Convert the coordinates to cylindrical coordinates.
        double r     = std::sqrt(x*x+y*y);
        double theta = std::atan2(y,x);

        // Compute cos and sin.
        double c     = std::cos(theta);
        double s     = std::sin(theta);

        std::array<double,6> cylField,cartField;
        cylField = field_rep.ComputeFieldInTime(t,r,theta,z,0);

        // Convert the cylindrical components to
        cartField[0] = c*cylField[0]-s*cylField[1];
        cartField[1] = s*cylField[0]+c*cylField[1];
        cartField[2] = cylField[2];
        cartField[3] = c*cylField[3]-s*cylField[4];
        cartField[4] = s*cylField[3]+c*cylField[4];
        cartField[5] = cylField[5];

        return cartField;
    }

protected:
    TemporalEMFieldMeshless<FieldRepresentation>  &  field_rep;
};

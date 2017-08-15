# ------------------------------- Information ------------------------------- #
# Author:       Joey Dumont                  <joey.dumont@gmail.com>          #
# Created:      2017-08-03                                                    #
# Description:  Creates initial conditions that are distributed with respect  #
#               to pair production probability.                               #
# --------------------------------------------------------------------------- #


# --------------------------- Modules Importation --------------------------- #
import sys as sys
import argparse as ap
from scipy.constants import codata
import numpy as np
from numpy import sin as sin
from numpy import cos as cos
from math import pow as pow
import xml.etree.ElementTree as ET
import time
# -- Import user modules.
sys.path.append("strattoanalysis/")
import GenerateInitialConditions
import Analysis3D

def main():
    """
    This class generates initial conditions for classical charged particle
    trajectory calculations that are distributed according to the probability
    of creating a pair at a given point in space-time.

    We proceed as follows:
        • Read the HDF5 file containing the temporal electromagnetic field.
        • Find the maximum value of the pair production integrand across the
            whole space and scale it.
        • Draw initial conditions for each time step by generating a random
            number r for each spatial point on the grid and instantiating a particle
            if scaledProbability > r until nParticles have been drawn.

    Author: Joey Dumont     <joey.dumont@gmail.com>
    """

    # Start the timer.
    start_time = time.perf_counter()

    # Command line arguments
    parser = ap.ArgumentParser(description="Generate the initial conditions for the simulation.")
    parser.add_argument("--directory", type=str, default=None,
                        help="Directory containing the HDF5 with the electromagnetic field data.")
    parser.add_argument("--fileTime",  type=str, default="Field_reflected_time.hdf5",
                        help="Name of the HDF5 file with the temporal electromagnetic field data.")
    parser.add_argument("--fileFreq",  type=str, default="Field_reflected.hdf5")
    parser.add_argument("--config",    type=str, default="configStrattoLinear.xml",
                        help="Name of the XML configuration file.")
    parser.add_argument("--radial",    type=bool,default=False,
                        help="Specifies the polarization of the beam.")

    # Parse arguments
    args = parser.parse_args()

    # Chosen shape for the simulation
    configFile = args.config

    # Parse arguments
    tree = ET.parse(configFile)
    config = tree.getroot()

    wavelength_element = config.find("./integration_salamin/lambda")
    if (wavelength_element == None):
        wavelength_element = config.find("./spectrum/lambda_c")
    wavelength = float(wavelength_element.text)

    # -- Electronic units.
    EL_UNITS_LENGTH = 2.0 * np.pi / wavelength
    EL_UNITS_TIME   = 2.0 * np.pi * codata.value('speed of light in vacuum') / wavelength

    numpart = int(config.find("./generate_initial_conditions/numpart").text)

    # We instantiate the Analysis3D object of interest and find the maximum pair density.
    if args.radial:
        pairProductionAnalysis = AnalysisRadial.AnalysisRadial(freq_field=args.directory+args.fileFreq, time_field=args.directory+args.fileTime)
    else:
        pairProductionAnalysis = Analysis3D.Analysis3D(freq_field=args.directory+args.fileFreq, time_field=args.directory+args.fileTime)

    maxIndices, maxDensities = pairProductionAnalysis.FindMaximumValues(pairProductionAnalysis.PairDensity)
    maxDensity               = np.amax(maxDensities)

    # We now scale the number of particles per time slices w.r.t the relative strength of maxDensities.
    maxDensities = maxDensities / np.sum(np.abs(maxDensities))
    numpart_slices = maxDensities * numpart

    # We prepare the arrays that will hold the initial positions and time values.
    x = np.empty((numpart))
    y = np.empty((numpart))
    z = np.empty((numpart))
    t = np.empty((numpart))

    for i in range(pairProductionAnalysis.size_time):
        particle_counter = 0
        loop_counter     = 0
        loop_max         = numpart
        while (particle_counter < numpart_slices[i] and loop_counter < loop_max):

            # -- Uniform numbers for this temporal slice.
            random_numbers = np.random.uniform(size=pairProductionAnalysis.size_flat)

            # -- Pair density for this slice.
            pairDensity = pairProductionAnalysis.PairDensityTime(i)/maxDensity
            print(pairDensity)

            for j in range(pairProductionAnalysis.size_flat):
                if (particle_counter >= numpart_slices[i]):
                    break
                if (pairDensity.flat[j] > random_numbers[j]):
                    indices = np.unravel_index(j, pairDensity.shape)

                    if args.radial:
                        r     = pairProductionAnalysis.coord_r[indices[0]]*pairProductionAnalysis.UNIT_LENGTH
                        theta = 2.0*np.pi*np.random.random()
                        z_si  = pairProductionAnalysis.coord_z[indices[1]]*pairProductionAnalysis.UNIT_LENGTH
                    else:
                        r       = pairProductionAnalysis.coord_r[indices[0]]*pairProductionAnalysis.UNIT_LENGTH
                        theta   = pairProductionAnalysis.coord_theta[indices[1]]
                        z_si    = pairProductionAnalysis.coord_z[indices[2]]*pairProductionAnalysis.UNIT_LENGTH

                    t_si = pairProductionAnalysis.time[i]*pairProductionAnalysis.UNIT_TIME

                    x[particle_counter] = r*np.cos(theta) * EL_UNITS_LENGTH
                    y[particle_counter] = r*np.sin(theta) * EL_UNITS_LENGTH
                    z[particle_counter] = z_si            * EL_UNITS_LENGTH
                    t[particle_counter] = t_si            * EL_UNITS_TIME
                    particle_counter    = particle_counter + 1

                    print("Allocated particle {} at x={}, y={}, z={} at t={}".format(particle_counter,x[particle_counter-1],y[particle_counter-1],z[particle_counter-1],t[particle_counter-1]))
                loop_counter += 1

    px = py = pz = 0.0

    # Create file
    of = open("init_conds.txt",'w')

    for pid in range(numpart): # Loop on particle indices
        # Write positions
        of.write(str(x[pid]) + " " + str(y[pid]) + " " + str(z[pid]) + " ")

        # Write momenta
        of.write(str(px) + " " + str(py) + " " + str(pz) + " ")

        # Write initial time.
        of.write(str(t[pid]) + "\n")

    # -- Stop timer.
    end_time = time.perf_counter()

    print(end_time-start_time)

if __name__ == "__main__":
    main()

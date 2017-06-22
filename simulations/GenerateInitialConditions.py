import sys as sys
import numpy as np
import argparse as ap
from numpy import sin as sin
from numpy import cos as cos
from math import pow as pow
import xml.etree.ElementTree as ET
#from IPython import embed

def placeParticlesInSphere(sphere_radius, wavelength, numpart):
    # Sphere radius in electronic units (converted from SI units, i.e. meters)
    R = 2.0 * np.pi * sphere_radius / wavelength

    # Generate uniformly distributed values of phi, cos theta and radius (in
    # spherical coordinates. radius is between 0 and 1)
    phi = np.random.uniform(low=0.0,high=2.0*np.pi,size=numpart)
    costheta = np.random.uniform(low=-1.0,high=1.0,size=numpart)
    radius = np.random.uniform(high=1.0,size=numpart)

    x = np.empty((numpart))
    y = np.empty((numpart))
    z = np.empty((numpart))
    for pid in range(numpart): # Loop on particle indices
        theta = np.arccos( costheta[pid] )
        r = R * pow( radius[pid], 1./3. ) # In electronic units (use R)

        x[pid] = r * sin( theta) * cos( phi[pid] )
        y[pid] = r * sin( theta) * sin( phi[pid] )
        z[pid] = r * cos( theta )

    return x, y, z


def placeParticlesInCylinder(base_radius, half_height, wavelength, numpart):
    # base radius in electronic units (converted from SI units, i.e. meters)
    B = 2.0 * np.pi * base_radius / wavelength
    # half height in electronic units
    H = 2.0 * np.pi * half_height / wavelength
    print("I am in cylinder")
    # Generate uniformly distributed values of cos theta, position in z and radius (in
    # cylindrical coordinates. radius is between 0 and 1)
    theta = np.random.uniform(low=0.0,high=2.0*np.pi,size=numpart)
    radius = np.random.uniform(high=1.0,size=numpart)
    randpos = np.random.uniform(low=-1.0,high=1.0,size=numpart)

    x = np.empty((numpart))
    y = np.empty((numpart))
    z = np.empty((numpart))
    for pid in range(numpart): # Loop on particle indices
        r = B * radius[pid]
        h = H * randpos[pid]

        x[pid] = r * cos(theta[pid])
        y[pid] = r * sin(theta[pid])
        z[pid] = 0.0

    return x, y, z

def placeParticlesInLine(initial_z, half_length, wavelength, numpart):
    # half length and initial position in z in electronic units
    L = 2.0 * np.pi * half_length / wavelength
    Z = 2.0 * np.pi * initial_z / wavelength

    # Generate uniformly distributed values of position in x
    randpos = np.random.uniform(low=-1.0,high=1.0,size=numpart)

    x = np.empty((numpart))
    y = np.empty((numpart))
    z = np.empty((numpart))
    for pid in range(numpart): # Loop on particle indices
        l = L * randpos[pid]

        x[pid] = l
        y[pid] = 0.0
        z[pid] = Z

    return x, y, z

def main():
    """
    Generate initial conditions for particle calculations.
    Takes input arguments in SI units on the command line
    and returns appropriate file with
    intial conditions in electronic units

    Author: D. Gagnon <denis.gagnon@emt.inrs.ca>
    """

    # Command line arguments
    parser = ap.ArgumentParser(description="Generate the initial conditions for the simulation.")
    parser.add_argument("--shape",  type=str, default="sphere",
                        help="Shape the initial conditions will be generated as.")

    # Parse arguments
    args = parser.parse_args()

    # Choosen shape for the simulation
    shape = args.shape

    # Parse arguments
    tree = ET.parse('config.xml')
    config = tree.getroot()
    wavelength = ""
    pz = ""
    numpart = ""
    sphere_radius = ""
    base_radius= ""
    half_height= ""
    initial_z= ""
    half_length= ""
    for child in config:
        if child.tag == "generate_initial_conditions":
            if  child[0].tag == "pz" \
            and child[1].tag == "numpart" \
            and child[2].tag == "sphere" \
            and child[3].tag == "cylinder" \
            and child[4].tag == "line":
                pz = float(child[0].text)
                numpart = int(child[1].text)
                if child[2][0].tag == "sphere_radius":
                    sphere_radius = float(child[2][0].text)
                else:
                    print("Wrong order of arguments in sphere tag")
                    sys.exit()
                if  child[3][0].tag == "base_radius" \
                and child[3][1].tag == "half_height":
                    base_radius = float(child[3][0].text)
                    half_height = float(child[3][1].text)
                else:
                    print("Wrong order of arguments in cylinder tag")
                    sys.exit()
                if  child[4][0].tag == "initial_z" \
                and child[4][1].tag == "half_length":
                    initial_z = float(child[4][0].text)
                    half_length = float(child[4][1].text)
                else:
                    print("Wrong order of arguments in line tag")
                    sys.exit()
            else:
                print("Wrong order of arguments in generate_initial_conditions tag")
                sys.exit()
        if child.tag == "integration_salamin":
            if  child[0].tag == "lambda":
                wavelength = float(child[0].text)
            else:
                print("Wrong order of arguments in integration_salamin tag")
                sys.exit()
    if wavelength == "" or pz == "" or numpart == "" \
    or sphere_radius == "" or base_radius == "" \
    or half_height == "" or initial_z == "" or half_length == "":
        print("Can't find generate_initial_conditions tag, integration_salamin tag or missing/empty argument")
        sys.exit()

    # Momentum values for px and py
    px = 0.0 
    py = 0.0
    # Momentum value for pz in electronic units (converted from eV/c)
    pz = pz / 0.5109989461e06 # Denominator is electron mass in eV/c^2

    if shape == "cylinder":
        x, y, z = placeParticlesInCylinder(base_radius, half_height, wavelength, numpart)
    elif shape == "line":
        x, y, z = placeParticlesInLine(initial_z, half_length, wavelength, numpart)
    else:
        x, y, z = placeParticlesInSphere(sphere_radius, wavelength, numpart)

    # Create file
    of = open("init_conds.txt",'w')

    for pid in range(numpart): # Loop on particle indices
        # Write positions
        of.write(str(x[pid]) + " " + str(y[pid]) + " " + str(z[pid]) + " ")

        # Write momenta
        of.write(str(px) + " " + str(py) + " " + str(pz) + '\n')

if __name__ == "__main__":
    main()

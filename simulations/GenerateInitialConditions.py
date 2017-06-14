import sys as sys
import numpy as np
from numpy import sin as sin
from numpy import cos as cos
from math import pow as pow
import xml.etree.ElementTree as ET
#from IPython import embed

def main():
    """
    Generate initial conditions for particle calculations.
    Takes input arguments in SI units on the command line
    and returns appropriate file with
    intial conditions in electronic units

    Author: D. Gagnon <denis.gagnon@emt.inrs.ca>
    """

    # Parse arguments
    tree = ET.parse('config.xml')
    config = tree.getroot()
    wavelength = ""
    pz = ""
    numpart = ""
    radius = ""
    for child in config:
        if child.tag == "generate_initial_conditions":
            if  child[0].tag == "pz" \
            and child[1].tag == "numpart" \
            and child[2].tag == "radius":
                pz = float(child[0].text)
                numpart = int(child[1].text)
                radius = float(child[2].text)
            else:
                print("Wrong order of arguments in generate_initial_conditions tag")
                sys.exit()
        if child.tag == "integration_salamin":
            if  child[0].tag == "lambda":
                wavelength = float(child[0].text)
            else:
                print("Wrong order of arguments in integration_salamin tag")
                sys.exit()
    if wavelength == "" or pz == "" or numpart == "" or radius == "":
        print("Can't find generate_initial_conditions tag, integration_salamin tag or missing argument")
        sys.exit()

    # Momentum values for px and py
    px = 0.0 
    py = 0.0
    # Momentum value for pz in electronic units (converted from eV/c)
    pz = pz / 0.5109989461e06 # Denominator is electron mass in eV/c^2

    # Sphere radius in electronic units (converted from SI units, i.e. meters)
    R = 2.0 * np.pi * radius / wavelength

    # Generate uniformly distributed values of phi, cos theta and radius (in
    # spherical coordinates. radius is between 0 and 1)
    phi = np.random.uniform(low=0.0,high=2.0*np.pi,size=numpart)
    costheta = np.random.uniform(low=-1.0,high=1.0,size=numpart)
    radius = np.random.uniform(high=1.0,size=numpart)

    # Create file
    of = open("init_conds.txt",'w')

    for pid in range(numpart): # Loop on particle indices
        theta = np.arccos( costheta[pid] )
        r = R * pow( radius[pid], 1./3. ) # In electronic units (use R)

        x = r * sin( theta) * cos( phi[pid] )
        y = r * sin( theta) * sin( phi[pid] )
        z = r * cos( theta )

        # Write positions
        of.write(str(x) + " " + str(y) + " " + str(z) + " ")

        # Write momenta
        of.write(str(px) + " " + str(py) + " " + str(pz) + '\n')

if __name__ == "__main__":
    main()

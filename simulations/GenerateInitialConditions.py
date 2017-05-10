import numpy as np
from numpy import sin as sin
from numpy import cos as cos
from math import pow as pow
import argparse
from IPython import embed

def main():
    """
    Generate initial conditions for particle calculations.
    Takes input arguments in SI units on the command line 
    and returns appropriate file with
    intial conditions in electronic units
    
    Author: D. Gagnon <denis.gagnon@emt.inrs.ca>
    """
    
    numpart = 25 # Number of particles
    R = 2.0 # Sphere radius
    
    # Command line arguments
    parser = argparse.ArgumentParser(description='Generate initial conditions for particle calculations.')
    parser.add_argument('--lambda',  type=float, default=1e06,
                        help='Laser wavelength in SI units (meters)')
    parser.add_argument('--pz', type=float,   default=0.0, 
                        help='Initial value of momentum in z in eV/c')
    parser.add_argument('--numpart', type=int,   default=10 
                        help='Number of initial conditions')
    parser.add_argument('--outfile', type=str,   default='init_conds.txt', 
                        help='Name of output file')

    # Parse arguments
    args = parser.parse_args()
   
    # Momentum values for px and py
    px = 0.0; py = 0.0; 
    # Momentum value for pz in electronic units (converted from eV/c)
    pz = args.pz
   
    # Generate random X positions
    phi = np.random.uniform(low=0.0,high=6.28,size=numpart)
    costheta = np.random.uniform(low=-1.0,high=1.0,size=numpart)
    radius = np.random.uniform(high=1.0,size=numpart)
    
    # Create file
    of = open('init_conds.txt','w')
    
    for pid in range(numpart):
        theta = np.arccos( costheta[pid] )
        r = R * pow( radius[pid], 1./3. )
        
        x = r * sin( theta) * cos( phi[pid] )
        y = r * sin( theta) * sin( phi[pid] )
        z = -6.0 + r * cos( theta )
        
        # Write momenta
        of.write(str(px) + " " + str(py) + " " + str(pz) + " ")
        
        # Write positions
        of.write(str(x) + " " + str(y) + " " + str(z) + '\n')
    
if __name__ == "__main__":
    main()

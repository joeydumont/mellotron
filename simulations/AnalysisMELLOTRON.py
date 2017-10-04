# ------------------------------- Information ------------------------------- #
# Author:       Joey Dumont                    <joey.dumont@gmail.com>        #
# Created       Oct 3rd, 2017                                                 #
# Description:  Collection of functions designed to analyze the collective    #
#               statistics of accelerated charges.                            #
# Dependencies: - NumPy                                                       #
#               - SciPy                                                       #
#               - H5Py                                                        #
#               - Matplotlib                                                  #
# --------------------------------------------------------------------------- #

# --------------------------- Modules Importation --------------------------- #
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
import scipy.signal as signal
import scipy.integrate as integration
import argparse
import h5py
import time

# --------------------------- Function Definition --------------------------- #




# ------------------------------- SCRATCH PAD ------------------------------- #

if __name__ == "__main__":

	# -- Command line arguments
	parser = ap.ArgumentParser(description="Produce the plots for a given .hdf5 output file."
	                    "If no file provided, the global.hdf5 file will be considered.")
	parser.add_argument("--directory", type=str,   default="./",
	                    help="Target directory containing a .hdf5 file.")
	parser.add_argument("--file", type=str,   default="global.hdf5",
	                    help=".hdf5 file with 1 particle to generate the trajectory plot with.")
	
	# -- Find the particle with the largest final gamma and plot its electric field.
	
	




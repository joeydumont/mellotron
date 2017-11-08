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
import argparse as ap
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import constants
import scipy.signal as signal
import scipy.integrate as integration
import argparse
import h5py as hp
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

	# -- Parse arguments
	args = parser.parse_args()

	# -- Make sure target directory wnds with a trailing slash.
	directory = args.directory
	if (not directory.endswith("/")):
		directory += "/"

	# -- Open the global hdf5 file.
	hdf5File = args.file
	#globalModel = findFileInDir(directory, hdf5File)
	globalModelFile = hp.File(directory + hdf5File, "r")

	# -- Plot the particles in the their initial position, and colour them according
	# -- to their final gamma values.
	initialX   = globalModelFile["/position"][0,:,0]
	initialY   = globalModelFile["/position"][0,:,1]
	initialZ   = globalModelFile["/position"][0,:,2]
	finalGamma = globalModelFile["/gamma"][-1,:]

	finalGamma = (4.00 * constants.physical_constants['atomic mass unit-electron volt relationship'][0]) * (finalGamma - 1.0)

	figInitialConditionsGamma = plt.figure(figsize=(8,4))
	figInitialConditionsGamma.subplots_adjust(hspace=0.3,wspace=0.3)

	axTranverse = plt.subplot2grid((1,2), (0,0))

	im = plt.scatter(initialX, initialY, c=finalGamma, norm=colors.LogNorm(),alpha=0.5)
	axTranverse.set_aspect('equal')
	axTranverse.set_xlabel('Initial $x$ position')
	axTranverse.set_ylabel('Initial $y$ position')

	divider = make_axes_locatable(axTranverse)
	cax = divider.append_axes("right", size="5%", pad=0.1)
	cbar = plt.colorbar(im, cax=cax)


	axLongitudinal = plt.subplot2grid((1,2), (0,1))
	im = plt.scatter(initialZ, initialY, c=finalGamma,norm=colors.LogNorm(), alpha=0.5)
	axLongitudinal.set_aspect('equal')
	axLongitudinal.set_xlabel('Initial $z$ position')

	divider = make_axes_locatable(axLongitudinal)
	cax = divider.append_axes("right", size="5%", pad=0.1)
	cbar = plt.colorbar(im, cax=cax)

	# -- Find the particle with the largest final gamma and plot its electric field.
	times       = globalModelFile["/times"][:]
	maxGammaIdx = np.argmax(globalModelFile["/gamma"][-1,:])
	maxChi      = globalModelFile["/chi"][:,maxGammaIdx]
	maxGammaTime= globalModelFile["/gamma"][:,maxGammaIdx]
	efield_Ex   = globalModelFile["electric_field"][:,maxGammaIdx,0]
	efield_Ey   = globalModelFile["electric_field"][:,maxGammaIdx,1]
	efield_Ez   = globalModelFile["electric_field"][:,maxGammaIdx,2]

	print(np.amax(efield_Ex))

	#print(efield_Ez)
	#print(efield_Ey)
	#print(efield_Ex)
	#print(maxChi)


	figElectricFieldMaxGamma = plt.figure(figsize=(4,3))
	ax1 = figElectricFieldMaxGamma.add_subplot(111)
	plt.plot(times, efield_Ex)
	plt.plot(times, efield_Ey)
	plt.plot(times, efield_Ez)
	#plt.plot(times, efield_Ex**2+efield_Ey**2+efield_Ez**2, 'k--')

	ax2 = ax1.twinx()
	plt.plot(times, (4.00 * constants.physical_constants['atomic mass unit-electron volt relationship'][0])*(maxGammaTime-1), 'k--')
	#plt.plot(times,maxChi)

	ax2.set_xlim((-150,150))

	figHistogramFinalGamma = plt.figure(figsize=(4,3))
	ax1 = figHistogramFinalGamma.add_subplot(111)

	finalGammaSorted = np.sort(finalGamma)

	print(finalGammaSorted)

	clipIndex = -1
	finalGammaHistogram = finalGammaSorted[0:clipIndex]

	plt.hist(finalGammaHistogram,
					 bins = 10 ** np.linspace(np.log10(finalGammaHistogram.min()), np.log10(finalGammaHistogram.max()), int(np.ceil(1.5*np.sqrt(finalGammaHistogram.size))) ))
	ax1.set_xscale('log')

	plt.show()


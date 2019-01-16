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
	parser.add_argument("--particlemass", type=float,   default=4.0, help="Ion mass in atomic mass units (u)")
	parser.add_argument("--lowenergycutoff", type=float,   default=0.0, help="Energy cutoff for direction plot (eV)")
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
	
	particleMass    = args.particlemass
	lowEnergyCutoff = args.lowenergycutoff

	# -- Plot the particles in the their initial position, and colour them according
	# -- to their final kinetic energy.
	initialX   = globalModelFile["/position"][0,:,0]
	initialY   = globalModelFile["/position"][0,:,1]
	initialZ   = globalModelFile["/position"][0,:,2]
	finalGamma = globalModelFile["/gamma"][-1,:]

	finalKinetic = (particleMass * constants.physical_constants['atomic mass unit-electron volt relationship'][0]) * (finalGamma - 1.0)

	figInitialConditionsKinetic = plt.figure(figsize=(8,4))
	figInitialConditionsKinetic.subplots_adjust(hspace=0.3,wspace=0.3)

	axTranverse = plt.subplot2grid((1,2), (0,0))

	im = plt.scatter(initialX, initialY, c=finalKinetic, norm=colors.LogNorm(),alpha=0.5)
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
	
	# -- Plot the final direction of particles, and colour them according
	# -- to their final kinetic energy.
	p_x_final         = globalModelFile["momentum"][-1,:,0]
	p_y_final         = globalModelFile["momentum"][-1,:,1]
	p_z_final         = globalModelFile["momentum"][-1,:,2]
	
	pr = np.sqrt(p_x_final*p_x_final + p_y_final*p_y_final + p_z_final*p_z_final)
	phi   = np.degrees(np.arctan2(p_y_final,p_x_final))
	theta = np.degrees(np.arccos(p_z_final/pr))

	figDirectionKinetic = plt.figure(figsize=(10,5))
	#figDirectionKinetic.subplots_adjust(hspace=0.3,wspace=0.3)

	axDirection = plt.gca()
	
	phiReduced   = []
	thetaReduced = []
	kinReduced   = []
	for i in range(len(p_x_final)):
		if finalKinetic[i] >= lowEnergyCutoff:
			phiReduced.append(phi[i])
			thetaReduced.append(theta[i])
			kinReduced.append(finalKinetic[i])
			
	im = plt.scatter(phiReduced, thetaReduced, c=kinReduced, norm=colors.LogNorm(),alpha=0.5, s=6)
			
	axDirection.set_aspect('equal')
	axDirection.set_xlabel('$\phi$ (degrees)')
	axDirection.set_ylabel('$\\theta$ (degrees)')
	axDirection.xaxis.set_ticks([-180, -90, 0, 90, 180])
	axDirection.yaxis.set_ticks([0, 45, 90, 135, 180])

	divider = make_axes_locatable(axDirection)
	cax = divider.append_axes("right", size="5%", pad=0.1)
	cbar = plt.colorbar(im, cax=cax)


	# -- Find the particle with the largest final kinetic energy and plot its electric field.
	times       = globalModelFile["/times"][:]
	maxGammaIdx = np.argmax(globalModelFile["/gamma"][-1,:])
	maxChi      = globalModelFile["/chi"][:,maxGammaIdx]
	maxGammaTime= globalModelFile["/gamma"][:,maxGammaIdx]
	efield_Ex   = globalModelFile["electric_field"][:,maxGammaIdx,0]
	efield_Ey   = globalModelFile["electric_field"][:,maxGammaIdx,1]
	efield_Ez   = globalModelFile["electric_field"][:,maxGammaIdx,2]
	p_x         = globalModelFile["momentum"][:,maxGammaIdx,0]
	p_y         = globalModelFile["momentum"][:,maxGammaIdx,1]
	p_z         = globalModelFile["momentum"][:,maxGammaIdx,2]

	# -- Compute force on particle.
	f_x = np.zeros_like(p_x)
	f_y = np.zeros_like(p_y)
	f_z = np.zeros_like(p_z)

	for i in range(1,p_x.size-1):
		f_x[i] = (p_x[i+1]-p_x[i-1])/(2*(times[1]-times[0]))
		f_y[i] = (p_y[i+1]-p_y[i-1])/(2*(times[1]-times[0]))
		f_z[i] = (p_z[i+1]-p_z[i-1])/(2*(times[1]-times[0]))

	mag_force = 2*np.pi * constants.c / 800e-9 * constants.m_e * constants.c * np.sqrt(f_x**2+f_y**2+f_z**2)

	plt.figure()
	plt.plot(times*800e-9 / (2*np.pi * constants.c),mag_force/(7344*constants.m_e))
	plt.gca().set_xlabel("Time")
	plt.gca().set_ylabel("Acceleration (m/s^2)")
	plt.savefig("ForceVSTime.pdf", bbox_inches='tight', dpi=500)


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
	plt.plot(times, (particleMass * constants.physical_constants['atomic mass unit-electron volt relationship'][0])*(maxGammaTime-1), 'k--')
	#plt.plot(times,maxChi)

	ax2.set_xlim((-150,150))

	figHistogramFinalKinetic = plt.figure(figsize=(4,3))
	ax1 = figHistogramFinalKinetic.add_subplot(111)

	finalKineticSorted = np.sort(finalKinetic)

	print(finalKineticSorted)

	clipIndex = -1
	finalKineticHistogram = finalKineticSorted[0:clipIndex]

	plt.hist(finalKineticHistogram,
					 bins = 10 ** np.linspace(np.log10(finalKineticHistogram.min()), np.log10(finalKineticHistogram.max()), int(np.ceil(1.5*np.sqrt(finalKineticHistogram.size))) ))
	ax1.set_xscale('log')
	ax1.set_xlabel('Kinetic energy (eV)')
	ax1.set_ylabel('Number of particles')
	

	plt.show()


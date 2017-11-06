# ------------------------------- Information ------------------------------- #
# Author:       Justine Pepin                  <justine.pepin@emt.inrs.ca>    #
# Created:      2017-06-06                                                    #
# Description:  Produces multiple plots to analyse MELLOTRON data.            #
# --------------------------------------------------------------------------- #

# --------------------------- Modules Importation --------------------------- #
import os as os
import argparse as ap
import h5py as hp
import sys as sys
from scipy import constants
import numpy.fft as fft
import scipy.fftpack as weird_fft
import scipy.integrate as integration
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# ------------------------------ Configuration ------------------------------ #
# Plot parameters
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'

# ---------------------------- Utility Functions ---------------------------- #
def findFileInDir(directory, fileName):
    """
    Determines if the specified HDF5 exists in the given directory.
    """
    # Find the .hdf5 file
    globalModel = ""
    for file in os.listdir(directory):
        if file == fileName:
            globalModel = file
            break

    if globalModel == "":
        print("It seems like the folder you gave doesn\'t have the specified .hdf5 file in it.")
        sys.exit()

    return globalModel

def setPositionsPlotLabels(graph):
    """
    Write the titles of the axis for the position plot.
    """
    graph.set_xlabel(r"X Axis")
    graph.set_ylabel(r"Y Axis")
    graph.set_zlabel(r"Z Axis")
    graph.set_title(r"Positions of Particles")

def createOneParticleTrajectory(directory, hdf5File, nTimeSteps, globalModelGroup):
    """
    Create a position plot (trajectories in space) in one-particle situations.

    This function ends the script because there is no need for polar plots in such a case.

    Trajectories are all different colors for better visualization purpose.
    """
    # Open position hdf5 table
    globalModelPositions = globalModelGroup["position"]
    # Initiate arrays of good size.
    xp = np.empty((nTimeSteps))
    yp = np.empty((nTimeSteps))
    zp = np.empty((nTimeSteps))
    # Stepping through the TimeSteps.
    for i in range(nTimeSteps):
        xp[i] = globalModelPositions[i, 0]
        yp[i] = globalModelPositions[i, 1]
        zp[i] = globalModelPositions[i, 2]

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot(xp, yp, zs=zp, lw=0.5)
    setPositionsPlotLabels(ax)

    plt.savefig(directory + os.path.splitext(hdf5File)[0] + ".pdf")
    sys.exit()

def createPositionsPlot(globalModelPositions, nParticles, nTimeSteps, directory):
    """
    Create a position plot (trajectories in space) in many-particles situations.
    Trajectories are all different colors for better visualization purpose.
    """
    # Initiate arrays of good size.
    xp = np.empty((nParticles, nTimeSteps))
    yp = np.empty((nParticles, nTimeSteps))
    zp = np.empty((nParticles, nTimeSteps))
    # Stepping through the Particles and the TimeSteps.
    for j in range(nParticles):
        for i in range(nTimeSteps):
            xp[j, i] = globalModelPositions[i, j, 0]
            yp[j, i] = globalModelPositions[i, j, 1]
            zp[j, i] = globalModelPositions[i, j, 2]

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    for j in range(nParticles):
        ax.plot(xp[j], yp[j], zs=zp[j], lw=0.5)
    setPositionsPlotLabels(ax)

    plt.savefig(directory + "positionsPlot.pdf")

def createPolarGammaPlot(globalModelMomentums, globalModelGamma, nParticles, nTimeSteps, directory, ionmode, ionmass, L):
    """
    Create a polar plot (gamma, times of flight and polar distributions) in many-particles situations.
    Divisions in histograms are all different colors for better visualization purpose.
    """
    # Initiate arrays of good size.
    r = np.empty((nParticles))
    gam = np.empty((nParticles))
    theta = np.empty((nParticles))
    phi = np.empty((nParticles))
    px_list = np.empty((nParticles))
    minGamma = 1.0
    maxGamma = 0.0
    isGammaSuperiorToOne = True

    # Define some constants if ion mode if turned on
    a = 1.0
    b = 0.0
    if ionmode:
        b = - 1.0
        a = ionmass * constants.physical_constants['atomic mass unit-electron volt relationship'][0]

    # Ion mass in electron mass units
    Mass = ionmass/constants.physical_constants['electron mass in u'][0]

    # Stepping through the Particles
    for j in range(nParticles):
            px = globalModelMomentums[nTimeSteps - 1, j, 0]
            py = globalModelMomentums[nTimeSteps - 1, j, 1]
            pz = globalModelMomentums[nTimeSteps - 1, j, 2]
            #gam[j] = np.sqrt(1 + (px**2 + py**2 + pz**2)/Mass**2 ) # Calculate gamma
            gam[j] = globalModelGamma[nTimeSteps - 1, j] # Take value in file
            if gam[j] <= minGamma:
                isGammaSuperiorToOne = False
            if r[j] > maxGamma:
                maxGamma = r[j]

            r[j] = a * (gam[j] + b)
            theta[j] = np.arctan2(px,pz)
            phi[j] = np.arctan2(px,py)

            px_list[j] = px

    f = plt.figure()
    ax1 = f.add_subplot(221, projection='polar')
    ax2 = f.add_subplot(222, projection='polar')
    ax3 = f.add_subplot(223)
    ax4 = f.add_subplot(224)
    ax1.grid(True)
    ax1.set_rmax(r.max())
    ax2.grid(True)
    ax2.set_rmax(r.max())
    ax3.grid(True)
    ax4.grid(True)
    ax1.scatter(theta, r)
    ax1.set_rlim(0)
    ax1.set_rscale('log')
    ax2.scatter(phi, r)
    ax2.set_rlim(0)
    ax2.set_rscale('log')
    print(r)

    if ionmode:
        ax1.set_ylabel(r'$E_k$ [eV]', labelpad=30)
        ax2.set_ylabel(r'$E_k$ [eV]', labelpad=30)
    else:
        ax1.set_ylabel(r'$\gamma$', labelpad=30)
        ax2.set_ylabel(r'$\gamma$', labelpad=30)

    ax1.set_title(r"zx plane ($\theta$)", va='bottom')
    ax2.set_title(r"yx plane ($\phi$)", va='bottom')


    # Histogram of gamma with colourful bands
    N, bins, patches = ax3.hist(r,  bins=int(np.ceil(1.5*np.sqrt(nParticles))), color=[0.8, 0.8, 0.2])
    for i in range(len(patches)):
        patches[i].set_facecolor((np.random.random(1)[0], np.random.random(1)[0], np.random.random(1)[0]))
    ax3.set_ylabel(r"Number of particles")


    # Time of flight histogram
    if ionmode:
        tx = np.zeros_like(gam[px_list>0])
        tx = (Mass*L*gam[px_list>0])/(px_list[px_list>0]*constants.c)
        ax3.set_xlabel(r'$E_k$ [eV]')
        N, bins, patches = ax4.hist(tx, bins = 10 ** np.linspace(np.log10(tx.min()), np.log10(tx.max()), int(np.ceil(1.5*np.sqrt(nParticles))) ), color=[0.8, 0.8, 0.2])
        for i in range(len(patches)):
            patches[i].set_facecolor((np.random.random(1)[0], np.random.random(1)[0], np.random.random(1)[0]))
        ax4.set_xscale('log')
    else:
        tx = np.zeros_like(gam[px_list<0])
        tx = (L*gam[px_list<0])/(px_list[px_list<0]*constants.c)
        ax3.set_xlabel(r'$\gamma$')
        N, bins, patches = ax4.hist(tx, bins =int(np.ceil(1.5*np.sqrt(nParticles))), color=[0.8, 0.8, 0.2])
        for i in range(len(patches)):
            patches[i].set_facecolor((np.random.random(1)[0], np.random.random(1)[0], np.random.random(1)[0]))

    ax4.set_xlabel(r"$t_x$ [s]")
    if len(tx) != 0:
        f.text(1, 0, "Fastest time of flight: %.5e s" % tx.min(), ha='right', fontdict=None)
    #f.text(0.5, 0, "Detector distance: %.5f m" % L, ha='center', fontdict=None)

    if ionmode:
        f.text(0, 0, "Particle mass is %.5f u" % ionmass, fontdict=None)

    plt.savefig(directory + "polarGammaPlots.pdf")

def main():
    """
    Produce the plots for a given .hdf5 output file.

    Author: Justine Pepin <justine.pepin@emt.inrs.ca>
    """

    # Command line arguments
    parser = ap.ArgumentParser(description="Produce the plots for a given .hdf5 output file."
                        "If no file provided, the global.hdf5 file will be considered.")
    parser.add_argument("--directory", type=str,   default="./",
                        help="Target directory containing a .hdf5 file.")
    parser.add_argument("--file", type=str,   default="global.hdf5",
                        help=".hdf5 file with 1 particle to generate the trajectory plot with.")
    parser.add_argument("--nTimeSteps", type=int,   default=0,
                        help="Number of timeSteps to plot on. If 0, the maximal number of timeSteps will be used. Default is the maximum.")
    parser.add_argument("--ion", type=bool,   default=False,
                        help="Switch for 'ion' mode in plots")
    parser.add_argument("--ionmass", type=float,   default=4.0,
                        help="Ion mass in atomic mass units (u)")
    parser.add_argument("--L", type=float,   default=0.03,
                        help="Distance of detector placed in x direction, in meters")
    parser.add_argument("--lw", type=bool, default=False,
                        help="Switch to produce plots related to the Liénard-Wiechert fields.")

    # Parse arguments
    args = parser.parse_args()

    # Target directory
    directory = args.directory
    if( not directory.endswith("/")):
        directory += "/"

    # Target file
    hdf5File = args.file
    if hdf5File != "global.hdf5":
        globalModel = findFileInDir(directory, hdf5File)
    else:
        globalModel = findFileInDir(directory, "global.hdf5")

    # Determine exactly how many TimeSteps there are.
    globalModelFile = hp.File(directory + globalModel, "r")
    globalModelGroup = globalModelFile.require_group("/")
    globalModelTimes = globalModelGroup["times"]
    nTimeSteps = globalModelTimes.len()

    # Target nTimeSteps
    targetTime = args.nTimeSteps
    if targetTime != 0 and targetTime < nTimeSteps:
        nTimeSteps = targetTime

    # If only one particle, no need to do rest of main
    if hdf5File != "global.hdf5":
        createOneParticleTrajectory(directory, hdf5File, nTimeSteps, globalModelGroup)

    # Determine exactly how many Particles there are.
    globalModelPositions = globalModelGroup["position"]
    nParticles = globalModelPositions.shape[1]

    # Create positions plot
    createPositionsPlot(globalModelPositions, nParticles, nTimeSteps, directory)

    # Create polar, gamma and time of flight plot
    globalModelGamma = globalModelGroup["gamma"]
    globalModelMomentums = globalModelGroup["momentum"]
    ionmode = args.ion
    ionmass = args.ionmass
    L = args.L
    createPolarGammaPlot(globalModelMomentums, globalModelGamma, nParticles, nTimeSteps, directory, ionmode, ionmass, L)

    # -- We produce the plots related to the Liénard-Wiechert fields.
    # -- We are interested in
    # --    • the distribution of power over the spectrum (dU/domega),
    # --    • the distribution of the number of photons over the spectrum (dN/domega)
    # --    • the distribution of energy over the sphere (dU/dOmega),
    # --    • the distribution of the number of photons over the sphere (dN/dOmega).
    # -- The rest will have to be taken care of via ParaView.
    if (args.lw):
        #-- First, we read the fields from the global HDF5 file.
        # -- /todo Add error catching.
        globalModelLWFields = globalModelGroup["lienard-wiechert-fields"]
        globalModelLWEField = globalModelLWFields["electric_field"]
        globalModelLWBField = globalModelLWFields["magnetic_field"]

        # -- The first step is to compute the FFT of the electric and
        # -- magnetic fields over each point of the sphere.
        padded_length  = weird_fft.next_fast_len(globalModelLWEField.shape[0])
        dt             = globalModelTimes[1]-globalModelTimes[0]
        fft_freqs      = 2*np.pi*fft.rfftfreq(padded_length, d=dt)
        fft_length     = fft_freqs.size

        lw_e_field_fft = np.zeros((fft_length,globalModelLWEField.shape[1],globalModelLWEField.shape[2],globalModelLWEField.shape[3]), dtype=complex)
        lw_b_field_fft = np.zeros((fft_length,globalModelLWEField.shape[1],globalModelLWEField.shape[2],globalModelLWEField.shape[3]), dtype=complex)

        for i in range(globalModelLWEField.shape[1]):
            for j in range(globalModelLWEField.shape[2]):
                for k in range(globalModelLWEField.shape[3]):
                    lw_e_field_fft[:,i,j,k] = dt*(fft.rfft(globalModelLWEField[:,i,j,k], n=padded_length))
                    lw_b_field_fft[:,i,j,k] = dt*(fft.rfft(globalModelLWBField[:,i,j,k], n=padded_length))

        # -- To compute dU/domega, we compute a surface integral over the observation sphere.
        # -- It is simpler to do this frequency per frequency.
        # -- We first compute the integrand.
        lw_field_power_spectrum_integrand = np.empty((lw_e_field_fft.shape[0],
                                                      lw_e_field_fft.shape[1],
                                                      lw_e_field_fft.shape[2]))

        # -- FOR TESTING PURPOSES ONLY
        lw_theta = np.linspace(0.0,   np.pi, lw_e_field_fft.shape[1])
        lw_phi   = np.linspace(0.0, 2*np.pi, lw_e_field_fft.shape[2])
        lw_radius= 1e8

        for i in range(lw_field_power_spectrum_integrand.shape[0]):
            for j in range(lw_field_power_spectrum_integrand.shape[1]):
                for k in range(lw_field_power_spectrum_integrand.shape[2]):
                    # -- For convenience, read the field components individually.
                    Ex = lw_e_field_fft[i,j,k,0]
                    Ey = lw_e_field_fft[i,j,k,1]
                    Ez = lw_e_field_fft[i,j,k,2]
                    Bx = np.conj(lw_b_field_fft[i,j,k,0])
                    By = np.conj(lw_b_field_fft[i,j,k,1])
                    Bz = np.conj(lw_b_field_fft[i,j,k,2])

                    integrand_x_comp = np.real(Ey*Bz-Ez*By)*np.sin(lw_theta[j])*np.cos(lw_phi[k])
                    integrand_y_comp = np.real(Ez*Bx-Ex*Bz)*np.sin(lw_theta[j])*np.sin(lw_phi[k])
                    integrand_z_comp = np.real(Ex*By-Ey*Bx)*np.cos(lw_theta[j])

                    lw_field_power_spectrum_integrand[i,j,k] = np.sin(lw_theta[j])*(integrand_x_comp+integrand_y_comp+integrand_z_comp)

        # -- We can now perform the angular integration.
        lw_field_power_spectrum = np.zeros((lw_e_field_fft.shape[0]))

        for i in range(lw_field_power_spectrum.shape[0]):
            lw_field_power_spectrum[i] = 4*np.pi*lw_radius**2*integration.simps(integration.simps(lw_field_power_spectrum_integrand[i], x=lw_phi), x=lw_theta)

        # -- Conversion to number of photons (requires proper units for frequencies).

        # -- To compute dU/dOmega, we perform an integral over the frequencies (also requires proper units).


        # -- PLOTTING
        plt.figure()
        axPowerSpectrum = plt.subplot2grid((2,2), (0,0))
        plt.semilogy(fft_freqs, lw_field_power_spectrum)
        axPowerSpectrum.set_xlabel("Frequency")
        axPowerSpectrum.set_ylabel("Spectral power")

        axTimeDependence = plt.subplot2grid((2,2), (0,1))
        plt.plot(globalModelTimes[:], globalModelLWEField[:,12,6,0])
        axTimeDependence.set_xlabel("Time")
        axTimeDependence.set_ylabel("$E_x$ along z axis")

        plt.savefig(directory + "Lienard-Wiechart-Fields.pdf", bbox_inches='tight')

if __name__ == "__main__":
    main()

import os as os
import argparse as ap
import h5py as hp
import sys as sys
from scipy import constants

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Plot parameters
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'


def findFileInDir(directory, fileName):
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
    graph.set_xlabel(r"X Axis")
    graph.set_ylabel(r"Y Axis")
    graph.set_zlabel(r"Z Axis")
    graph.set_title(r"Positions of Particles")

def createOneParticleTrajectory(directory, hdf5File, nTimeSteps, globalModelGroup):
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
    
    plt.savefig(directory + os.path.splitext(hdf5File)[0] + ".eps")
    sys.exit()

def createPositionsPlot(globalModelPositions, nParticles, nTimeSteps, directory):
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

    plt.savefig(directory + "positionsPlot.eps")

def createPolarGammaPlot(globalModelMomentums, globalModelGamma, nParticles, nTimeSteps, directory, ionmode, ionmass, L):
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
    ax1.set_rmax(maxGamma)
    ax2.grid(True)
    ax2.set_rmax(maxGamma)
    ax3.grid(True)
    ax4.grid(True)
    ax1.scatter(theta, r)
    ax2.scatter(phi, r)

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
    tx = np.zeros_like(gam[px_list>0])
    if ionmode:
        tx = (Mass*L*gam[px_list>0])/(px_list[px_list>0]*constants.c) 
        ax3.set_xlabel(r'$E_k$ [eV]')
    else:
        tx = (L*gam[px_list>0])/(px_list[px_list>0]*constants.c) 
        ax3.set_xlabel(r'$\gamma$')
    
    N, bins, patches = ax4.hist(tx, bins = 10 ** np.linspace(np.log10(tx.min()), np.log10(tx.max()), int(np.ceil(1.5*np.sqrt(nParticles))) ), color=[0.8, 0.8, 0.2])
    for i in range(len(patches)):
        patches[i].set_facecolor((np.random.random(1)[0], np.random.random(1)[0], np.random.random(1)[0]))

    ax4.set_xscale('log')
    ax4.set_xlabel(r"$t_x$ [s]")
    f.text(1, 0, "Fastest time of flight: %.5e s" % tx.min(), ha='right', fontdict=None)
    #f.text(0.5, 0, "Detector distance: %.5f m" % L, ha='center', fontdict=None)
    
    if ionmode:
        f.text(0, 0, "Particle mass is %.5f u" % ionmass, fontdict=None)
        
    plt.savefig(directory + "polarGammaPlots.eps")

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

    # Create polar chi plot
    globalModelGamma = globalModelGroup["gamma"]
    globalModelMomentums = globalModelGroup["momentum"]
    ionmode = args.ion
    ionmass = args.ionmass
    L = args.L
    createPolarGammaPlot(globalModelMomentums, globalModelGamma, nParticles, nTimeSteps, directory, ionmode, ionmass, L)
    

if __name__ == "__main__":
    main()

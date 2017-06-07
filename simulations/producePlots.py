import os as os
import argparse as ap
import h5py as hp
import sys as sys

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Plot parameters
mpl.rcParams['ps.useafm'] = True
mpl.rcParams['pdf.use14corefonts'] = True
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'

def createPositionsPlot(globalModelPositions, nParticles, nTimeSteps, directory):
    # Initiate arrays of good size.
    xp = np.empty((nParticles, nTimeSteps))
    yp = np.empty((nParticles, nTimeSteps))
    zp = np.empty((nParticles, nTimeSteps))


    # Stepping through the Particles and the TimeSteps.
    for j in range(nParticles):
        for i in range(nTimeSteps):
            xp[j, i] = globalModelPositions[i, 0, j]
            yp[j, i] = globalModelPositions[i, 1, j]
            zp[j, i] = globalModelPositions[i, 2, j]

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    for j in range(nParticles):
        ax.plot(xp[j], yp[j], zs=zp[j], lw=0.5)
    ax.set_xlabel(r"X Axis")
    ax.set_ylabel(r"Y Axis")
    ax.set_zlabel(r"Z Axis")
    ax.set_title(r"Positions of Particles")

    plt.savefig(directory + "positionsPlot.eps")


def createPolarGammaPlot(globalModelMomentums, globalModelGamma, nParticles, nTimeSteps, directory):
    r = np.empty((nParticles))
    theta = np.empty((nParticles))
    phi = np.empty((nParticles))
    maxGamma = 0.0

    for j in range(nParticles):
            r[j] = globalModelGamma[nTimeSteps - 1, 0, j]
            if r[j] > maxGamma:
                maxGamma = r[j]
            px = globalModelMomentums[nTimeSteps - 1, 0, j]
            py = globalModelMomentums[nTimeSteps - 1, 1, j]
            pz = globalModelMomentums[nTimeSteps - 1, 2, j]
            theta[j] = np.arctan2(px,pz)
            phi[j] = np.arctan2(px,py)

    f = plt.figure()
    ax1 = f.add_subplot(221, projection='polar')
    ax2 = f.add_subplot(222, projection='polar')
    ax3 = f.add_subplot(223)
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax1.scatter(theta, r)
    ax2.scatter(phi, r)

    ax1.set_title(r"zx plane ($\theta$)", va='bottom')
    ax1.set_ylabel(r"gamma", labelpad=30)
    ax2.set_title(r"yx plane ($\phi$)", va='bottom')
    ax2.set_ylabel(r"gamma", labelpad=30)

    N, bins, patches = ax3.hist(r, bins=int(np.ceil(1.5*np.sqrt(nParticles))), color=[0.8, 0.8, 0.2])
    for i in range(len(patches)):
        patches[i].set_facecolor((np.random.random(1)[0], np.random.random(1)[0], np.random.random(1)[0]))

    ax3.set_xlabel(r"gamma")
    ax3.set_ylabel(r"Number of particles")
    plt.savefig(directory + "polarGammaPlots.eps")

def main():
    """
    Produce the plots for a given global.hdf5 file.

    Author: Justine Pepin <justine.pepin@emt.inrs.ca>
    """

    # Command line arguments
    parser = ap.ArgumentParser(description="Produce the plots for a given global.hdf5 file.")
    parser.add_argument("--directory", type=str,   default="./",
                        help="Target directory containing global.hdf5 file")

    # Parse arguments
    args = parser.parse_args()

    # Target directory
    directory = args.directory
    if( not directory.endswith("/")):
        directory += "/"

    # Find the global.hdf5 file
    globalModel = ""
    for file in os.listdir(directory):
        if file == "global.hdf5":
            globalModel = file
            break

    if globalModel == "":
        print("It seems like the folder you gave doesn\'t have a global.hdf5 file in it.")
        sys.exit()
    
    # Determine exactly how many TimeSteps there are.
    globalModelFile = hp.File(directory + globalModel, "r")
    globalModelGroup = globalModelFile.require_group(globalModel)
    globalModelTimes = globalModelGroup["times"]
    nTimeSteps = globalModelTimes.len()

    # Determine exactly how many Particles there are.
    globalModelPositions = globalModelGroup["position"]
    nParticles = globalModelPositions.shape[2]

    # Create positions plot
    createPositionsPlot(globalModelPositions, nParticles, nTimeSteps, directory)

    # Create polar chi plot
    globalModelGamma = globalModelGroup["gamma"]
    globalModelMomentums = globalModelGroup["momentum"]
    createPolarGammaPlot(globalModelMomentums, globalModelGamma, nParticles, nTimeSteps, directory)
    

if __name__ == "__main__":
    main()

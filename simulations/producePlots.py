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
            xp[j, i] = globalModelPositions[i, 0, j]
            yp[j, i] = globalModelPositions[i, 1, j]
            zp[j, i] = globalModelPositions[i, 2, j]

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    for j in range(nParticles):
        ax.plot(xp[j], yp[j], zs=zp[j], lw=0.5)
    setPositionsPlotLabels(ax)

    plt.savefig(directory + "positionsPlot.eps")

def calculateCOnV(r, nParticles):
    for j in range (nParticles):
        if r[j] == 1.0:
            print("Division by 0... exiting.")
            sys.exit()
        r[j] = 1.0/np.sqrt(1.0 - 1.0/np.power(r[j], 2))
    return r

def createPolarGammaPlot(globalModelMomentums, globalModelGamma, nParticles, nTimeSteps, directory, a, b):
    # Initiate arrays of good size.
    r = np.empty((nParticles))
    gam = np.empty((nParticles))
    theta = np.empty((nParticles))
    phi = np.empty((nParticles))
    minGamma = 1.0
    maxGamma = 0.0
    isGammaSuperiorToOne = True
    # Stepping through the Particles.
    for j in range(nParticles):
            r[j] = a * (globalModelGamma[nTimeSteps - 1, 0, j] + b)
            gam[j] = globalModelGamma[nTimeSteps - 1, 0, j]
            if gam[j] <= minGamma:
                isGammaSuperiorToOne = False
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
    ax4 = f.add_subplot(224)
    ax1.grid(True)
    ax1.set_rmax(maxGamma)
    ax2.grid(True)
    ax2.set_rmax(maxGamma)
    ax3.grid(True)
    ax4.grid(True)
    ax1.scatter(theta, r)
    ax2.scatter(phi, r)
    gammaLabel = "*(gamma + "
    gammaLabel = "(" + str(a) + ")" + gammaLabel + "(" + str(b) + "))"
    ax1.set_title(r"zx plane ($\theta$)", va='bottom')
    ax1.set_ylabel(gammaLabel, labelpad=30)
    ax2.set_title(r"yx plane ($\phi$)", va='bottom')
    ax2.set_ylabel(gammaLabel, labelpad=30)

    # Histogram of gamma with colourful bands
    N, bins, patches = ax3.hist(r, bins=int(np.ceil(1.5*np.sqrt(nParticles))), color=[0.8, 0.8, 0.2])
    for i in range(len(patches)):
        patches[i].set_facecolor((np.random.random(1)[0], np.random.random(1)[0], np.random.random(1)[0]))

    # Histogram of c/v with colourful bands
    if isGammaSuperiorToOne:
        gam = calculateCOnV(gam, nParticles)
        N, bins, patches = ax4.hist(gam, bins=int(np.ceil(1.5*np.sqrt(nParticles))), color=[0.8, 0.8, 0.2])
        for i in range(len(patches)):
            patches[i].set_facecolor((np.random.random(1)[0], np.random.random(1)[0], np.random.random(1)[0]))

    ax3.set_xlabel(gammaLabel)
    ax3.set_ylabel(r"Number of particles")
    ax4.set_xlabel(r"c/v")
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
    parser.add_argument("--a", type=float,   default=1.0,
                        help="a value in a*(gamma + b).")
    parser.add_argument("--b", type=float,   default=0.0,
                        help="b value in a*(gamma + b).")

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
    globalModelGroup = globalModelFile.require_group(globalModel)
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
    nParticles = globalModelPositions.shape[2]

    # Create positions plot
    createPositionsPlot(globalModelPositions, nParticles, nTimeSteps, directory)

    # Create polar chi plot
    globalModelGamma = globalModelGroup["gamma"]
    globalModelMomentums = globalModelGroup["momentum"]
    a = args.a
    b = args.b
    createPolarGammaPlot(globalModelMomentums, globalModelGamma, nParticles, nTimeSteps, directory, a, b)
    

if __name__ == "__main__":
    main()

import os as os
import argparse as ap
import h5py as hp
import sys as sys

def writeAttribute1(of, n, nTimeSteps, nParticles,nameOfAttribute):
    """
    Writes into the global.xdmf file an attribute that only has one dimension such as chi or gamma.
    """
    of.write('<Attribute Name="' + str(nameOfAttribute) + '" Center="Node">\n')
    of.write('<DataItem Format="HDF" ItemType="HyperSlab" Dimensions="1 ' + str(nParticles) + '">\n')
    of.write('<DataItem Dimensions="3 2" NumberType="Int">\n')
    of.write(str(n) + ' 0\n')
    of.write('1 1\n')
    of.write('1 ' + str(nParticles) + '\n')
    of.write('</DataItem>\n')
    of.write('<DataItem Name="' + str(nameOfAttribute) + '" Format="HDF" NumberType="Float" Dimensions="' + str(nTimeSteps) + ' ' + str(nParticles) + '">\n')
    of.write('global.hdf5:/' + str(nameOfAttribute) + '\n')
    of.write('</DataItem>\n')
    of.write('</DataItem>\n')
    of.write('</Attribute>\n')

def writeAttribute3(of, n, nTimeSteps, nParticles, nameOfAttribute):
    """
    Writes into the global.xdmf file an attribute that only has one dimension such as the electric field or the magnetic field.
    """
    of.write('<Attribute Name="' + str(nameOfAttribute) + '" Center="Node">\n')
    of.write('<DataItem Format="HDF" ItemType="HyperSlab" Dimensions="1 ' + str(nParticles) + ' 3">\n')
    of.write('<DataItem Dimensions="3 3" NumberType="Int">\n')
    of.write(str(n) + ' 0 0\n')
    of.write('1 1 1\n')
    of.write('1 ' + str(nParticles) + ' 3\n')
    of.write('</DataItem>\n')
    of.write('<DataItem Name="' + str(nameOfAttribute) + '" Format="HDF" NumberType="Float" Dimensions="' + str(nTimeSteps) + ' ' + str(nParticles) + ' 3">\n')
    of.write('global.hdf5:/' + str(nameOfAttribute) + '\n')
    of.write('</DataItem>\n')
    of.write('</DataItem>\n')
    of.write('</Attribute>\n')

def generateXMF(directory, nParticles, nTimeSteps):
    """
    Creates the global.xdmf file and writes the timesteps list in it.
    For each of the timesteps, this function writes the position in 3D and the attributes of the particle.
    """
    # Initialize xdmf file
    of = open(directory + 'global.xdmf','w')
    of.write('<?xml version="1.0" ?>\n')
    of.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
    of.write('<Xdmf Version="2.1">\n')
    of.write('<Domain>\n')
    of.write('<Grid Name="Temporal Collection" GridType="Collection" CollectionType="Temporal">\n')

    # Write the timesteps list
    of.write('<Time TimeType="List">\n')
    of.write('<DataItem ItemType="HyperSlab" Dimensions="' + str(nTimeSteps) + '">\n')
    of.write('<DataItem Dimensions="3 1">\n')
    of.write('0\n')
    of.write('1\n')
    of.write(str(nTimeSteps) + '\n')
    of.write('</DataItem>\n')
    of.write('<DataItem Name="times" Format="HDF" NumberType="Float" Dimensions="' + str(nTimeSteps) + '">\n')
    of.write('global.hdf5:/times\n')
    of.write('</DataItem>\n')
    of.write('</DataItem>\n')
    of.write('</Time>\n')

    # For each timestep
    for n in range(0, nTimeSteps):
        # Declare a grid of points
        of.write('<Grid Name="timestep' + str(n) + '" GridType="Uniform">\n')
        of.write('<Topology TopologyType="Polyvertex" NodesPerElement="1" NumberOfElements="' + str(nParticles) + '"/>\n')
        of.write('<Geometry GeometryType="XYZ">\n')
        of.write('<DataItem ItemType="HyperSlab" Dimensions="1 ' + str(nParticles) + ' 3">\n')
        of.write('<DataItem Dimensions="3 3" NumberType="Int">\n')
        of.write(str(n) + ' 0 0\n')
        of.write('1 1 1\n')
        of.write('1 ' + str(nParticles) + ' 3\n')
        of.write('</DataItem>\n')
        of.write('<DataItem Name="position" Format="HDF" NumberType="Float" Dimensions="' + str(nTimeSteps) + ' ' + str(nParticles) + ' 3">\n')
        of.write('global.hdf5:/position\n')
        of.write('</DataItem>\n')
        of.write('</DataItem>\n')
        of.write('</Geometry>\n')

        # Write attributes
        # -- chi
        writeAttribute1(of, n, nTimeSteps, nParticles, "chi")
        # -- gamma
        writeAttribute1(of, n, nTimeSteps, nParticles, "gamma")
        # -- magnetic_field
        writeAttribute3(of, n, nTimeSteps, nParticles, "magnetic_field")
        # -- electric_field
        writeAttribute3(of, n, nTimeSteps, nParticles, "electric_field")
        # -- momentum
        writeAttribute3(of, n, nTimeSteps, nParticles, "momentum")

        # Close grid after writing attributes
        of.write('</Grid>\n')

    # Close all nametags to finalize and close file
    of.write('</Grid>\n')
    of.write('</Domain>\n')
    of.write('</Xdmf>\n')
    of.close()


def addToGlobal(directory, partialfile, globalGroup, nTimeSteps, n):
    """
    Reads the chi, gamma, magnetic, electric, momentum and position values
    from the single particle HDF5 output and accumulates the data in the
    global file.
    """
    # Fill the particle dimension in global hdf5 file
    partialFile = hp.File(directory + partialfile, "r")
    partialGroup = partialFile.require_group(partialfile)
    # -- chi
    partialData = partialGroup["chi"]
    globalData = globalGroup["chi"]
    globalData[0:nTimeSteps, n] = partialData[0:nTimeSteps]
    # -- gamma
    partialData = partialGroup["gamma"]
    globalData = globalGroup["gamma"]
    globalData[0:nTimeSteps, n] = partialData[0:nTimeSteps]
    # -- magnetic_field
    partialData = partialGroup["magnetic_field"]
    globalData = globalGroup["magnetic_field"]
    globalData[0:nTimeSteps, n, 0:3] = partialData[0:nTimeSteps, 0:3]
    # -- electric_field
    partialData = partialGroup["electric_field"]
    globalData = globalGroup["electric_field"]
    globalData[0:nTimeSteps, n, 0:3] = partialData[0:nTimeSteps, 0:3]
    # -- momentum
    partialData = partialGroup["momentum"]
    globalData = globalGroup["momentum"]
    globalData[0:nTimeSteps, n, 0:3] = partialData[0:nTimeSteps, 0:3]
    # -- position
    partialData = partialGroup["position"]
    globalData = globalGroup["position"]
    globalData[0:nTimeSteps, n, 0:3] = partialData[0:nTimeSteps, 0:3]

    # -- close partial hdf5 file
    partialFile.close()

def accumulateLWInGlobal(directory, partialfile, globalGroup, nTimeSteps):
    """
    Reads the Liénard-Wiechert fields from the single particle HDF5 output
    and accumulates the data in the global file.
    """
    # -- Open the single-particle file (group is the same as the filename).
    partialFile  = hp.File(directory + partialfile, "r")
    partialGroup = partialFile.require_group(partialfile+"/lienard-wiechert-fields")

    # -- Accumulate the fields.
    globalGroup["electric_field"][:] += partialGroup["electric_field"][:]
    globalGroup["magnetic_field"][:] += partialGroup["magnetic_field"][:]

    # -- Close the resources.
    partialFile.close()

def main():
    """
    Manage the .hdf5 output files to create a global .hdf5 file complemented
    with a .xdmf file to be read in Paraview.

    Author: Justine Pepin <justine.pepin@emt.inrs.ca>
    """

    # Command line arguments
    parser = ap.ArgumentParser(description="Manage the .hdf5 output files.")
    parser.add_argument("--directory", type=str,   default="./",
                        help="Target directory containing output files")
    parser.add_argument("--lw", type=bool, default=False,
                        help="Switch whether to sum the Liénard-Wiechert fields.")

    # Parse arguments
    args = parser.parse_args()

    # Target directory
    directory = args.directory
    if( not directory.endswith("/")):
        directory += "/"

    # Find a times model and calculate the number of particles
    timesModel = ""
    nParticles = 0
    for file in os.listdir(directory):
        if file.endswith(".hdf5"):
            if file != "global.hdf5":
                nParticles += 1
                timesModel = file

    if nParticles == 0 or timesModel == "":
        print("It seems like the folder you gave doesn\'t have hdf5 files in it.")
        sys.exit()

    # Determine exactly how many time steps there are.
    timesModelFile = hp.File(directory + timesModel, "r")
    timesModelGroup = timesModelFile.require_group(timesModel)
    timesModelTimes = timesModelGroup["times"]
    nTimeSteps = timesModelTimes.len()

    # Create canvas of global hdf5 file
    globalFile = hp.File(directory + "global.hdf5", "w")
    globalGroup = globalFile.require_group("/")

    # -- times
    globalGroup.copy(timesModelTimes, "times", "times", False, False, False, False, False)
    # -- chi
    globalGroup.create_dataset("chi", (nTimeSteps, nParticles), dtype="f8")
    # -- gamma
    globalGroup.create_dataset("gamma", (nTimeSteps, nParticles), dtype="f8")
    # -- magnetic_field
    globalGroup.create_dataset("magnetic_field", (nTimeSteps, nParticles, 3), dtype="f8")
    # -- electric_field
    globalGroup.create_dataset("electric_field", (nTimeSteps, nParticles, 3), dtype="f8")
    # -- momentum
    globalGroup.create_dataset("momentum", (nTimeSteps, nParticles, 3), dtype="f8")
    # -- position
    globalGroup.create_dataset("position", (nTimeSteps, nParticles, 3), dtype="f8")

    # Find all .hdf5 files in given directory
    n = -1
    for file in os.listdir(directory):
        if file.endswith(".hdf5"):
            if file != "global.hdf5":
                if n < nParticles:
                    n = n + 1
                else:
                    break
                addToGlobal(directory, file, globalGroup, nTimeSteps, n)

    # -- If LW is set, we create a lienard-wiechert-fields group in global.hdf5,
    # -- copy the theta and phi datasets from a single partial HDF5 file,
    # -- then accumulate the fields emitted by all particles.
    if (args.lw):
        # -- We copy the data to groups in the global hdf5 file.
        LW_ModelGroup   = timesModelGroup["lienard-wiechert-fields"]
        globalLWGroup = globalGroup.create_group("lienard-wiechert-fields")
        LW_ModelGroup.copy("phi",   globalLWGroup)
        LW_ModelGroup.copy("theta", globalLWGroup)

        # -- We now create the field datasets in the global file.
        globalLWGroup.create_dataset("electric_field", (nTimeSteps, LW_ModelGroup["electric_field"].shape[1], LW_ModelGroup["electric_field"].shape[2], 3), dtype=float, fillvalue=0.0)
        globalLWGroup.create_dataset("magnetic_field", (nTimeSteps, LW_ModelGroup["electric_field"].shape[1], LW_ModelGroup["electric_field"].shape[2], 3), dtype=float, fillvalue=0.0)

        n = -1
        for file in os.listdir(directory):
            if file.endswith(".hdf5"):
                if file != "global.hdf5":
                    if n < nParticles:
                        n = n + 1
                    else:
                        break
                    accumulateLWInGlobal(directory, file, globalLWGroup, nTimeSteps)


    # -- Close remaining resources.
    timesModelFile.close()
    globalFile.close()

    generateXMF(directory, nParticles, nTimeSteps)

if __name__ == "__main__":
    main()
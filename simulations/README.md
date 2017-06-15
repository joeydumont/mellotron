PROJECT DESCRIPTION:
MELLOTRON is a C++ template, header-only library aimed at the efficient
computation of charged particles trajectories in electromagnetic fields. Its main
design feature is the use of arbitrary electromagnetic field models through the use of
templates. It is based upon the Boost integration library, Boost.Numeric.Odeint.
Particles trajectories can be computed for any charge/mass ratio.



STARTING THE MELLOTRON:
1) Insure the .o are present in the simulations/ directory. 

2) Write your simulation parameters in the config.xml file without changing the file's layout.

3) Create a new directory that will contain your simulation results :

        mkdir <dirname>

   If you are missing inspiration for the name, you can always use the current date and time like that :

        mkdir $(date +%F--%T)

4) Copy the config.xml file into your new directory without changing the name of the file :

        cp config.xml <dirname>/

4.a) OPTIONAL add pre-calculated initial conditions and/or normalization constant to your new directory.
     They will be considered by the MELLOTRON and used instead of re-calculated. They must be named either

        init_conds.txt

     for the initial conditions or

        normalization_constant.txt

     for the normalization constant.

        cp init_conds.txt <dirname>/
        cp normalization_constant.txt <dirname>/

5) Launch the MELLOTRON :

        ./automaticStart.sh <dirname>/

5.a) OPTIONAL generate the plots for your simulation.

        python producePlots.py --directory <dirname>/



SCRIPTS API:
1) automaticStart.sh                            Launches the MELLOTRON if provided at least a dirname containing a config.xml file.
        usage : ./automaticStart.sh <dirname>/
        note : you might have to do a small chmod +x automaticStart.sh before to obtain permission to use the script.

2) ComputeNormalizationConstantSalaminLinear.o  Computes the normalization constant the MELLOTRON needs.
        usage : There is a config.xml in the same directory containing the right values of lambda, w0 and L.
                ./ComputeNormalizationConstantSalaminLinear.o
        note : don't forget to compile the most recent version of .cpp with compile.sh script.

3) GenerateInitialConditions.py                 Generates the x,y,z,px,py and pz of particles that will be simulated in the MELLOTRON
        usage : There is a config.xml in the same directory containing the right values of lambda, pz numpart and radius.
                python GenerateInitialConditions.py

4) IntegrationSalamin.o                         Calculates the trajectory of each initial condition
        usage : There is a config.xml in the same directory containing the right values of mass, Q, lambda, w0, L, energy, t_init, dt and nsteps.
                There is a init_conds.txt file at disposition.
                There is a normalization_constant.txt file in the same directory.
                cat path/to/init_conds.txt | parallel -j <number of jobs> --colsep " " ./IntegrationSalamin.o --init_conds {1} {2} {3} {4} {5} {6}

5) manageOutputs.py                             Generates the global.hdf5 file containing the data of all trajectories calculated in the MELLOTRON and its global.xdmf file.
        usage : python manageOutputs.py --nParticles <number of .hdf5 output files> --directory <dirname>/
                paraview <dirname>/global.xdmf

6) producePlots.py                              Generates the plots of trajectory and/or polar position and gamma after a given .hdf5 file.
        usage : python producePlots.py --directory <dirname>/ --file <name.hdf5> --nTimeSteps <int>
                default value for --file option is global.hdf5.
                default value for nTimeSteps is 0, which means all timeSteps will be used to generate the trajectory.
                epstopdf <dirname>/<NameOfGeneratedPlot>.eps <dirname>/<NameOfGeneratedPlot>.pdf
                xdg-open <dirname>/<NameOfGeneratedPlot>.pdf
        note : If given a global.hdf5 file, with python producePlots.py --directory <dirname>/, it will produce a positions plot with all the trajectories on.
               It will also produce a polar positions and gamma plot.
               If given a single particle .hdf5 file (e.g. 14595196478285746701.hdf5), it will only produce a positions plot with the single trajectory on.



CONTRIBUTORS:
Joey Dumont <joey.dumont@emt.inrs.ca>
Denis Gagnon <denis.gagnon@emt.inrs.ca>
Justine Pepin <justine.pepin@emt.inrs.ca>

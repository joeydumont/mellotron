PROJECT DESCRIPTION:
MELLOTRON is a C++ template, header-only library aimed at the efficient
computation of charged particles trajectories in electromagnetic fields. Its main
design feature is the use of arbitrary electromagnetic field models through the use of
templates. It is based upon the Boost integration library, Boost.Numeric.Odeint.
Particles trajectories can be computed for any charge/mass ratio.



STARTING THE MELLOTRON FROM A LOCAL COMPUTER :
1) Insure the .o are present in the simulations/ directory. 

2) Write your simulation parameters in the configSalamin.xml/configStrattoLinear.xml file without changing the file's layout.

3) Create a new directory that will contain your simulation results :

        mkdir <dirname>

   If you are missing inspiration for the name, you can always use the current date and time like that :

        mkdir $(date +%F--%T)

4) Copy the configuration file into your new directory without changing the name of the file :

        cp configSalamin.xml <dirname>/

4.a) OPTIONAL add pre-calculated initial conditions and/or normalization constant to your new directory.
     They will be considered by the MELLOTRON and used instead of re-calculated. They must be named either

        init_conds.txt

     for the initial conditions or

        normalization_constant.txt

     for the normalization constant.

        cp init_conds.txt <dirname>/
        cp normalization_constant.txt <dirname>/

5) Launch the MELLOTRON :

        ./automaticStart.sh {options | -j <int> | -s <cylinder, sphere, line> | -d <dirname> | -c <configSalamin.xml, configStrattoLinear.xml> }



STARTING THE MELLOTRON FROM MAMMOUTH :
1) Move to a directory onto the scratch. Insure you copied on the scratch a configuration file (either configSalamin.xml or configStrattoLinear.xml) and the mammouthStart.sh file.
        cd /mnt/parallel_scratch_mp2_wipe_on_december_2017/maclean/maclean_group/mellotron

2) Write your simulation parameters in the configSalamin.xml/configStrattoLinear.xml file without changing the file's layout. 
   Write the arguments in the mammouthStart file, plus change the PBS values for the submission. PBS node must equals to NNODES.

3) Create a new directory that will contain your simulation results :

        mkdir <dirname>

   If you are missing inspiration for the name, you can always use the current date and time like that :

        mkdir $(date +%F--%T)

4) Copy the configuration file into your new directory without changing the name of the file, and the mammouthStart.sh as well :

        cp configSalamin.xml <dirname>/
        cp mammouthStart.sh <dirname>/

4.a) OPTIONAL add pre-calculated initial conditions and/or normalization constant to your new directory.
     They will be considered by the MELLOTRON and used instead of re-calculated. They must be named either

        init_conds.txt

     for the initial conditions or

        normalization_constant.txt

     for the normalization constant.

        cp init_conds.txt <dirname>/

5) Launch the MELLOTRON :

        qsub mammouthStart.sh



SCRIPTS API:
1.a) automaticStart.sh                          Launches the MELLOTRON if provided at least a dirname containing a configSalamin.xml file.
        usage : ./automaticStart.sh -d <dirname>/ -j <number of jobs> -s <shape> -c <config.xml>
                default value for -d <dirname> is current directory.
                default value for -j <number of jobs> is two. 
                default value for -s <shape> is sphere. 
                default value for -c <config.xml> is configSalamin.xml.
        note :  you might have to do a small chmod +x automaticStart.sh before to obtain permission to use the script.
                The <shape> parameter takes a string that can be chosen between sphere, cylinder or line. If given something else, sphere will be used.
                The <config.xml> parameter takes a string that can be either be configSalamin.xml or configStrattoLinear.xml.

1.b) mammouthStart.sh                           Launches the MELLOTRON when on mammouth.
        usage : qsub mammouthStart.sh
        #PBS -l walltime=120:00:00 This line takes the maximum estimated time of execution in format hours:minutes:seconds.
        #PBS -l nodes=16:ppn=1     This line takes the number of nodes that will be used by mammouth for the simulation.
        note :  The NNODES parameter takes the same int than passed to the PBS nodes argument.
                The SHAPE parameter takes a string that can be chosen between sphere, cylinder or line.
                The CONFIG parameter takes a string that can be either be configSalamin.xml or configStrattoLinear.xml.

2) ComputeNormalizationConstantSalaminLinear.o  Computes the normalization constant the MELLOTRON needs.
        usage : There is a configSalamin.xml in the same directory containing the right values of lambda, w0 and L.
                ./ComputeNormalizationConstantSalaminLinear.o
        note :  don't forget to compile the most recent version of .cpp with compile.sh script or cmake . in parent directory called mellotron/.

3) GenerateInitialConditions.py                 Generates the x,y,z,px,py and pz of particles that will be simulated in the MELLOTRON
        usage : There is a configuration file in the same directory containing the right values of lambda/lambda_c, pz and numpart.
                The values for the shapes are also there : initial_z, half_length, half_height, base_radius and sphere_radius.
                python GenerateInitialConditions.py --shape <shape> --config <config.xml>
                default value for <shape> is sphere. Also, cylinder and line are valid strings that can be passed.
                default value for <config.xml> is configSalamin.xml. Also, configStrattoLinear.xml is a valid value.

4) IntegrationSalamin.o                         Calculates the trajectory of each initial condition
        usage : There is a configSalamin.xml in the same directory containing the right values of mass, Q, lambda, w0, L, energy, t_init, dt and nsteps.
                There is a init_conds.txt file at disposition.
                There is a normalization_constant.txt file in the same directory.

4) IntegrationStrattoLinear.o                   Calculates the trajectory of each initial condition
        usage : There is a configStrattoLinear.xml in the same directory.
                There is a init_conds.txt file at disposition.

5) manageOutputs.py                             Generates the global.hdf5 file containing the data of all trajectories calculated in the MELLOTRON and its global.xdmf file.
        usage : python manageOutputs.py --nParticles <number of .hdf5 output files> --directory <dirname>/
                paraview <dirname>/global.xdmf

6) producePlots.py                              Generates the plots of trajectory and/or polar position and gamma after a given .hdf5 file.
        usage : python producePlots.py --directory <dirname>/ --file <name.hdf5> --nTimeSteps <int> --ion <bool> --ionmass <float> --L <float>
                default value for --file option is global.hdf5.
                default value for --nTimeSteps option is 0, which means all timeSteps will be used to generate the trajectory.
                default value for --ion option is False
                default value for --ionmass option is 4.0.
                default value for --L option is 0.03.
                epstopdf <dirname>/<NameOfGeneratedPlot>.eps <dirname>/<NameOfGeneratedPlot>.pdf
                xdg-open <dirname>/<NameOfGeneratedPlot>.pdf
        note :  If given a global.hdf5 file, with python producePlots.py --directory <dirname>/, it will produce a positions plot with all the trajectories on.
                It will also produce a polar positions and gamma plot.
                If the `ion` flag is set to False, the gamma factor will be plotted.
                If the `ion` flag is set to True, the kinetic energy will be plotted instead of the gamma factor. This kinetic energy value will depend
                on the `ionmass` value, which should be in atomic mass units (u). The value of `L` is used in the time of flight plot. It corresponds to the
                distance (in meters) between an imaginary detector placed parallel to the x-axis and the origin of the coordinate system.



CONTRIBUTORS:
Joey Dumont <joey.dumont@emt.inrs.ca>
Denis Gagnon <denis.gagnon@emt.inrs.ca>
Justine Pepin <justine.pepin@emt.inrs.ca>

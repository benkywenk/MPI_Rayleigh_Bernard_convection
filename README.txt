2.5D and 3D Rayleigh Bernard convection within rotating frame. 

This code can be used so that simulations can be run on Carpathia, there are 4 important files in the repository:

    1: rayleigh_bernard.py, this is the relevant program to run if you want to perform a 2.5D rotating simulation.
       Within the file is more detailed instructions of how it works.

    2: 3D_rayleigh_bernard.py, this is similar to the first however unsurprisingly in full 3D. Again the file itself
       is commented to explain what each thing is doing.

    3: run_param_file.py, this is the file to change if you want to adjust the variables of the simulation. 
    
    4: merge.py, after a simulation has been completed this file should be run to merge the hdf5 files into one file.
       This process works with data produced from multiple cores or one core.

First step once you have ssh'd into Carpathia, is to type the two following commands:

    source /opt/Miniconda/miniconda37/etc/profile.d/conda.sh

    conda activate dedalus

Since dedalus is allready installed on Carpathia, this first command points to where the python libraries are and 
then opens the dedalus library. At this point, you can edit and run simulations however you want. As explained earlier,
the file you will edit the most is the parameters file. Within the ssh, you can use the following command to quickly
edit the parameter file from within the terminal:

    nano run_param_file.py 

The same can be done for the other files as well if you replace the file name with the corresponding file. I would 
advise that when you are writing longer pieces of code that you use the ssh functionality of vscode. When you do run
a simulation, you can run it using the following command:

    python3 rayleigh_benard.py

However this is only really done for a really low resolution simulation. The more useful command to use to run a 
simulation using mpi (on for example, 4 processors) can be done by typing:

    mpiexec -n 4 python3 rayleigh_benard.py

This command can be changed to select a higher or lower number of cores. Carpathia has 64 cores and realistically, 
for a 2.5D simulation 4 or 8 cores can be used. Don't set this much hugher in a 2.5D simulation as it is unnecessary.
For 3D simulations, 16 cores can be used however be warned that the time it takes to calculate timesteps increases alot
faster in 3D. on 16 cores, I would advise a resolution of either 64x64x32 or 64x64x64, which take approx 2.5-3 hrs and 10-12 hrs.
Once we have a working angular momentum spectral decomposition analysis tool then we can start to ramp up the resolution
and take a bit more time to make a bigger dataset (However realistically we need to just work on one dataset between us,
or mayve we would each make different simulations with a few different parameters to investigate how these things change the simulation).

BEFORE RUNNING:
make sure you have created a file structure that corresponds to how the data will be saved. In my case I created a folder called sim_data
and made my save_direc variable be equal to that. 

Once a simulation has been completed, the following command can be done:

    python3 merge.py

This doesn't need to be done using mpiexec because then each core would produce an indentical file, so it is not needed.

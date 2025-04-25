This repository provides various example LAMMPS simulation files as well as all scripts written by me as part of the project.


The 1_lammps_file_setup fodler contains scripts that help setup files for LAMMPS simulations:

    - script.py - A python script for the setup of a LAMMPS input .data file that is required to run a 
    simulation.
    - z1_lammps_script_generation.ipynb - A jupyter notebook that is able to easily write and create a LAMMPS
    input .in file that is required to run a simulation. There is also functionality to run a quick local 
    LAMMPS simulation here if desired.


The 2_lammps_simulation_inputs folder contains example input files for LAMMPS simulations:

    - input.in - An example LAMMPS input .in file.
    - job.sh - An example job submission .sh file.
    - lammps_input_data.data - An example LAMMPS input .data file for a CC3 cage in dichloromethane solvent.


The 3_analysis_scripts folder contains scripts related to the analysis of the simulations:

    - z2_simulation_stability_script.ipynb - A jupyter notebook that extracts simulation data from the LAMMPS
    output logfile and plots graphs against time.
    - z3_convert_lammstring_to_xyz.ipynb - A jupyter notebook that converts the output LAMMPS .lammstrj file
    to various .xyz files.
    - z4_pore_volumes_script.ipynb - A jupyter notebook that writes various .xyz files for visualisation, 
    while also producing analytic data for the simulation, including: 
        - pywindow pore volume
        - pywindow pore diameter
        - pywindow pore windows
        - convex hull volumes
    - z5_solvent_inside_or_outside_cage.ipynb - A jupyter notebook that writes various .xyz files for 
    visualisation, while also producing analytic data for the simulation, including:
        - the number of solvents inside the cage based on pywindow and convex hull definitions
        - solvent hopping rates based on pywindow and convex hull definitions.
    - z6_symmetry_script.ipynb - A jupyter notebook that writes various .xyz files for visualisation, while
    also producing analytic data for the simulation, including:
        - standard deviations of the angles and edge lengths of the underlying polyhedron of the cage
        - tritopic-ditopic-tritopic angles
        - dihedral angles
        - other code and ideas for alternative symmetry based analysis
    

The 4_conformer_energy_analysis folder contains files related to the analysis of the conformer energy of the 
cage over a series of frames of a simulation:

    - input_frame_0.in - An example LAMMPS input .in file for a single frame that is evaluated for cage 
    conformer energy analysis. A short minimisation is run here as it provides the initial energy of the 
    system (the cage), which is the conformer energy of the cage for a given frame.
    - SPE_job.sh - An example of job submission .sh file that loops through every evaluated frame for cage 
    conformer energy analysis.
    - SPE_script.py - A python script for the setup of the LAMMPS input .data files that are required to run a 
    simulation. One .data file is produced for each frame that is evaluated.
    - z8_SPE_script.ipynb - A jupyter notebook that extracts the conformer energies of the cage from the ouput
    LAMMPS logfiles and plots them against time.


The requirements.txt file contains the packages and the versions used for the .ipynb and .py files presented 
in this repository.


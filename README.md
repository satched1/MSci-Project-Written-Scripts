Within the master_scripts_file, this repository provides various example LAMMPS simulation files as well as all scripts written by me as part of the project.


he 1_lammps_file_setup folder contains scripts that help setup files for LAMMPS simulations:
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
    to various .xyz files, detailed below.
    - z4_pore_volumes_script.ipynb - A jupyter notebook that writes various .xyz files for visualisation, 
    while also producing analytic data for the simulation, including: 
        - pywindow pore volume
        - pywindow pore diameter
        - pywindow pore windows
        - convex hull volumes
    - z5_solvent_inside_or_outside_cage.ipynb - A jupyter notebook that writes various .xyz files, detailed 
    below, for visualisation, while also producing analytic data for the simulation, including:
        - the number of solvents inside the cage based on pywindow and convex hull definitions
        - solvent hopping rates based on pywindow and convex hull definitions.
    - z6_symmetry_script.ipynb - A jupyter notebook that writes various .xyz files, detailed below, for 
    visualisation, while also producing analytic data for the simulation, including:
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
    - z8_SPE_script.ipynb - A jupyter notebook that extracts the conformer energies of the cage from the output
    LAMMPS logfiles and plots them against time.


The requirements.txt file contains the packages and the versions used for the .ipynb and .py files presented 
in this repository.


Here are details on the .xyz files that can be outputted from the written scripts:
    - Files from the z3_convert_lammstring_to_xyz.ipynb:
        - data_all_atoms_unwrapped.xyz - xyz file containing coordinate information for all atoms for each 
        frame in an unwrapped coordinate style.
        - data_all_atoms_wrapped.xyz - xyz file containing coordinate information for all atoms for each 
        frame in a wrapped coordinate style.
        - data_cage_only_unwrapped.xyz - xyz file containing coordinate information for all cage atoms for 
        each frame in an unwrapped coordinate style.
        - data_cage_only_wrapped.xyz - xyz file containing coordinate information for all cage atoms for each 
        frame in a wrapped coordinate style.
        - data_solvent_only_unwrapped.xyz - xyz file containing coordinate information for all solvent atoms 
        for each frame in an unwrapped coordinate style.
        - data_solvent_only_wrapped.xyz - xyz file containing coordinate information for all solvent atoms 
        for each frame in a wrapped coordinate style.
    - Files from the z4_pore_volumes_script.ipynb:
        - coms_tritopic_bbs.xyz - xyz file containing centre of mass information for the tritopic building 
        blocks of the cage for each frame.
    - Files from the z5_solvent_inside_or_outside_cage.ipynb:
        - coms_tritopic_bbs.xyz - xyz file containing centre of mass information for the tritopic building 
        blocks of the cage for each frame.
        - Cu_all_solvents.xyz - xyz file containing the cage and all solvent molecules for each frame, 
        where each solvent molecule is represented by a Cu atom at its centre of mass.
        - Cu_ch_solvents.xyz - xyz file containing the cage and all solvent molecules that are defined as 
        being inside the cage for each frame based on a convex hull pore volume definition, where each 
        solvent molecule is represented by a Cu atom at its centre of mass.
        - Cu_large_radius_solvents.xyz - xyz file containing the cage and all solvent molecules that are 
        defined as being inside the cage for each frame based on a large threshold radius definition, where 
        each solvent molecule is represented by a Cu atom at its centre of mass.
        - Cu_pw_solvents.xyz - xyz file containing the cage and all solvent molecules that are defined as 
        being inside the cage for each frame based on a pywindow pore volume definition, where each solvent 
        molecule is represented by a Cu atom at its centre of mass.
    - Files from the z6_symmetry_script.ipynb:
        - coms_all_bbs.xyz - xyz file containing centre of mass information for all building blocks of the 
        cage for each frame.
        - coms_ditopic_bbs.xyz - xyz file containing centre of mass information for the ditopic building 
        blocks of the cage for each frame.
        - coms_tritopic_bbs.xyz - xyz file containing centre of mass information for the tritopic building 
        blocks of the cage for each frame.


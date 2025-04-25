# import libraries
from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np
import openmm

from openff.toolkit.topology import Molecule as molk
from openff.toolkit.topology import Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.units import unit
from openff.interchange import Interchange

import stk
from rdkit import Chem 
import time
from espaloma_charge import charge
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
from openff.interchange.components._packmol import UNIT_CUBE, pack_box
from openff.toolkit import  Molecule

import os

print("done")

openmm.Platform.getNumPlatforms()

toolkit_registry = EspalomaChargeToolkitWrapper() # function to help assign partial charges

def chiral_stereo_check(mol):
    '''
    Function to help stereochemistry assignments.
    '''
    Chem.AddHs(mol)
    Chem.SanitizeMol(mol)
    Chem.Kekulize(mol)
    Chem.FindMolChiralCenters(mol,includeUnassigned=True, force=True)
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True) #includeUnassigned
    Chem.AssignAtomChiralTagsFromStructure(mol)    
    return mol

def get_cage_frame_data(file_name, frame_number):
    '''
    Outputs the position matrix of the cage for a given frame.

    Inputs:
    file_name (str): path to the cage only unwrapped file
    frame_number (int): the frame to be extracted
    Returns:
    frame_data (list): list[cage_atom][coorindates] containing cage position information
    '''
    # open file and setup counts and conditions
    datafile = open(file_name, 'r')
    line_count = 0
    frame_count = 0 # count of lines after frame has been identified
    frame_id = False
    frame_data = [] # empty list to store cage frame information

    # loop through each line
    for line in datafile:

        # define number of atoms in cage on first line of file
        if line_count == 0:
            frame_atoms_range = range(0,int(line))
        line_count = line_count + 1

        # when condition is true, extract frame coordinates
        if frame_id == True:
            if frame_count in frame_atoms_range:
                line = line.split()
                frame_data.extend([[float(x) for x in line[1:]]])
            frame_count = frame_count + 1
        
        # set condition to true when the desired frame is encountered
        if line == f"frame {frame_number}\n":
            frame_id = True
    return frame_data


# loop through desired frames from the simulation to output LAMMPS input .data files for
for i in range(0,201,5):
    
    # build cage with stk
    bb1 = stk.BuildingBlock('CCC1=C(C(=C(C(=C1CN)CC)CN)CC)CN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock.init_from_file('dit_6.mol',functional_groups=[stk.AldehydeFactory()],)
    cage1 = stk.ConstructedMolecule(
                topology_graph=stk.cage.EightPlusTwelve(building_blocks={
                        bb2: range(8, 20),bb1: range(0, 8),},optimizer =stk.MCHammer(),),)
    

    # extract cage atomic coordinates and add 150 to each coordinate to place cage in the centre of a massive box
    centering = np.array(get_cage_frame_data('data_cage_only_unwrapped.xyz',i)) +150
    print(centering)
    
    # update the stk position matrix with the new coordinates
    cage1 = cage1.with_position_matrix(centering)


    topology1=Topology() # create openff topology object
    mol_list=[]      # set empty list for cage molecule     

    # cubic_box sets the simulation box size. Ensure this is comfortably larger than the size of the cage itself guarantee boundary effects when calculating energies are not encountered
    cubic_box = unit.Quantity(300 * np.eye(3), unit.angstrom)


    polymer = stk.BuildingBlock.to_rdkit_mol(cage1) # make rdkit molecule from stk
    molecule= chiral_stereo_check(polymer) # check rdkit stereochemistry


    cc = molk.from_rdkit(molecule,allow_undefined_stereo=False) # make openff molecule from rdkit molecule
    molk.assign_partial_charges(cc, partial_charge_method='espaloma-am1bcc', toolkit_registry=toolkit_registry) # assign partial charges with openff, using espaloma partial charge method
    topology1.add_molecule(cc) # add openff cage to topology object
    mol_list.append(cc) # add open cage to cage molecule list


    forcefield= ForceField("openff-2.1.0.offxml") # load openff forcefield
    print(mol_list)
    interchange = Interchange.from_smirnoff(forcefield,topology1,box=cubic_box,charge_from_molecules=mol_list) # FROM OPEN-FF TO OPEN-MM , box=cubic_box charge_from_molecules=[cc1,cc]
    system=interchange.to_openmm(combine_nonbonded_forces = True) # create openmm topology object



    # here is an optional openmm simulation, which can be included as a possible measure against the simulation initially blowing up. Not run here for conformer energy analysis.
    #5x small steps current state num steps trj freq data freq
    num_steps = 400  # number of integration steps to run
    trj_freq = 500  # number of steps per written trajectory frame
    data_freq = 400  # number of steps per written simulation statistics
    time_step = 1 * openmm.unit.femtoseconds  # simulation timestep normal 1 
    temperature = 298 * openmm.unit.kelvin  # simulation temperature normal 700
    friction = 1 / openmm.unit.picosecond  # friction constant   normal 1
    integrator = openmm.LangevinMiddleIntegrator(temperature, friction, time_step)
    simulation = interchange.to_openmm_simulation(integrator=integrator) # FROM OPEN-FF TO OPEN-MM
    #simulation.minimizeEnergy(tolerance=3, maxIterations=50) #new added !!! in some test the simulation atoms blew up so an #initial minimization can ammend this issue



    # here are various reporters to help output desired file and data types
    print("done")
    pdb_reporter = openmm.app.PDBReporter("trajectory.pdb", trj_freq) #,enforcePeriodicBox=True
    state_data_reporter = openmm.app.StateDataReporter(
        "data.csv",
        data_freq,
        step=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,
    )
    simulation.reporters.append(pdb_reporter) # reporter to report .pdb data
    simulation.reporters.append(state_data_reporter) # reporter to report data
    print("next")

    interchange.to_lammps(f"lammps_input_data_frame_{i}.data") # output lammps .data file

    os.rename(f"lammps_input_data_frame_{i}.data.lmp", f"lammps_input_data_frame_{i}.data") # rename .data file file as it is initially created with .lmp extension
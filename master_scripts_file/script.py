# import libraries
from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np
import openmm

from openff.toolkit.topology import Molecule as molk
from openff.toolkit.topology import Topology
from openff.toolkit.typing.engines.smirnoff import ForceField # smirnoff ff api part
from openff.units import unit
from openff.interchange import Interchange

import stk
from rdkit import Chem 
import time
import stk
import numpy as np
from espaloma_charge import charge
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
from openff.interchange.components._packmol import UNIT_CUBE, pack_box
from openff.toolkit import  Molecule

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


# build cage with stk
bb2 = stk.BuildingBlock('C1=C(C=C(C=C1C=O)C=O)C=O', [stk.AldehydeFactory()])
bb1 = stk.BuildingBlock(
    smiles='N[C@@H]1CCCC[C@H]1N',
    functional_groups=[stk.PrimaryAminoFactory()],
)
cage1 = stk.ConstructedMolecule(
            topology_graph=stk.cage.EightPlusTwelve(
                building_blocks={
                    bb1: range(8, 20),
                    bb2: range(0, 8),
                },
                optimizer =stk.MCHammer(),
                ),
        )


topology1=Topology() # create openff topology object
mol_list=[]    # set empty list for cage molecule     

molecules = [Molecule.from_smiles(smi) for smi in ["C(Cl)Cl"]]  #you can add additional solvent molecules here by appending to the list

# cubic_box sets the simulation box size. Ensure this is comfortably larger than the size of the cage itself
cubic_box = unit.Quantity(40 * np.eye(3), unit.angstrom)


# cage_temp is a displacement of the cage to centre of simulation box - we want to avoid any possible periodic boundary effects on setup
cage_temp= cage1.with_displacement([20,20,20])
polymer = stk.BuildingBlock.to_rdkit_mol(cage_temp) # make rdkit molecule from stk
molecule= chiral_stereo_check(polymer) # check rdkit stereochemistry

cc = molk.from_rdkit(molecule,allow_undefined_stereo=False) # make openff molecule from rdkit molecule
molk.assign_partial_charges(cc, partial_charge_method='espaloma-am1bcc', toolkit_registry=toolkit_registry) # assign partial charges with openff, using espaloma partial charge method
topology1.add_molecule(cc) # add openff cage to topology object
mol_list.append(cc) # add open cage to cage molecule list

# pack simulation box using openff interchange, packing the box using packmol
topology = pack_box(
    molecules=molecules,
    solute=topology1,
    number_of_copies=[798],  # number of solvent molecules
    target_density=1330 * unit.kilogram / unit.meter**3, # upper limit of density; note this line doesn't seem to do much since we already define simulation box volume and number of solvents
    box_shape=UNIT_CUBE,
)    


forcefield= ForceField("openff-2.1.0.offxml") # load openff forcefield
print(mol_list)
interchange = Interchange.from_smirnoff(forcefield,topology,box=cubic_box,charge_from_molecules=mol_list) # FROM OPEN-FF TO OPEN-MM , box=cubic_box charge_from_molecules=[cc1,cc]
system=interchange.to_openmm(combine_nonbonded_forces = True) # create openmm topology object



# here is an optional openmm simulation, which can be included as a possible measure against the simulation initially blowing up
#5x small steps current state num steps trj freq data freq
num_steps = 400  # number of integration steps to run
trj_freq = 500  # number of steps per written trajectory frame
data_freq = 400  # number of steps per written simulation statistics
time_step = 1 * openmm.unit.femtoseconds  # simulation timestep normal 1 
temperature = 298 * openmm.unit.kelvin  # simulation temperature normal 700
friction = 1 / openmm.unit.picosecond  # friction constant normal 1
integrator = openmm.LangevinMiddleIntegrator(temperature, friction, time_step) # define openmm langevin integrator
simulation = interchange.to_openmm_simulation(integrator=integrator) # FROM OPEN-FF TO OPEN-MM
#simulation.minimizeEnergy(tolerance=3, maxIterations=50)
print("done")



# here are various reporters to help output desired file and data types
pdb_reporter = openmm.app.PDBReporter("trajectory.pdb", trj_freq)
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
simulation.reporters.append(pdb_reporter) # reporter to report .pdb data, unecessary since we dont need the trajectory.pdb file?
simulation.reporters.append(state_data_reporter) # reporter to report data, uneceessary since we dont need the data.csv file?

print("next")

interchange.to_lammps("lammps_input_data.data") # output lammps .data file

import os
os.rename("lammps_input_data.data.lmp", "lammps_input_data.data") # rename .data file file as it is initially created with .lmp extension
# Plain MD simulation of benzene bound to T4 lysozyme

import pathlib

import gufe
from openfe import (
    ChemicalSystem,
    ProteinComponent,
    SmallMoleculeComponent,
    SolventComponent,
)
from openfe.protocols.openmm_md.plain_md_methods import PlainMDProtocol
from openff.units import unit

# Define the Chemical System

ligand = SmallMoleculeComponent.from_sdf_file("inputs/benzene.sdf")
solvent = SolventComponent(ion_concentration=0.15 * unit.molar)
protein = ProteinComponent.from_pdb_file("inputs/t4_lysozyme.pdb", name="t4-lysozyme")  # ty: ignore[invalid-argument-type]

system = ChemicalSystem(
    {"ligand": ligand, "protein": protein, "solvent": solvent},
    name=f"{ligand.name}_{protein.name}",
)

# Adjust Simulation Settings

settings = PlainMDProtocol.default_settings()
settings.simulation_settings.equilibration_length_nvt = 0.01 * unit.nanosecond
settings.simulation_settings.equilibration_length = 0.01 * unit.nanosecond
settings.simulation_settings.production_length = 0.02 * unit.nanosecond
settings.output_settings.checkpoint_interval = 0.02 * unit.nanosecond
settings.engine_settings.compute_platform = "CPU"
settings.solvation_settings.solvent_padding = 1.0 * unit.nanometer

print(settings)

# Create a Protocol

protocol = PlainMDProtocol(settings=settings)

# Running the Simulation

workdir = pathlib.Path("./")
dag = protocol.create(stateA=system, stateB=system, mapping=None)
dag_result = gufe.protocols.execute_DAG(
    dag,
    shared_basedir=workdir,
    scratch_basedir=workdir,
    keep_shared=True,
    n_retries=3,
)

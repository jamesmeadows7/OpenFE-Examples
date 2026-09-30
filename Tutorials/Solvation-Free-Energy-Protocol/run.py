# Absolute hydration free energy of benzene in water

import json
import pathlib

import gufe.tokenization
import openfe
from gufe.protocols import execute_DAG
from openfe.protocols.openmm_afe import AbsoluteSolvationProtocol
from openff.units import unit
from rdkit import Chem

# Load the Ligand

supp = Chem.SDMolSupplier("inputs/benzene.sdf", removeHs=False)
ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in supp]

# Create the Chemical Systems

solvent = openfe.SolventComponent()

systemA = openfe.ChemicalSystem(
    {  # ligand interacting with solvent
        "ligand": ligands[0],
        "solvent": solvent,
    },
    name=ligands[0].name,
)

systemB = openfe.ChemicalSystem({"solvent": solvent})  # ligand decoupled from solvent

# Creating the Protocol

settings = AbsoluteSolvationProtocol.default_settings()
settings.protocol_repeats = 1
settings.lambda_settings.lambda_elec = [
    0.0,
    0.26,
    0.5,
    0.75,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
]
settings.solvent_simulation_settings.equilibration_length = 10 * unit.picosecond
settings.solvent_simulation_settings.production_length = 500 * unit.picosecond
settings.solvent_engine_settings.compute_platform = "CPU"
settings.vacuum_simulation_settings.equilibration_length = 10 * unit.picosecond
settings.vacuum_simulation_settings.production_length = 500 * unit.picosecond
settings.vacuum_engine_settings.compute_platform = "CPU"

protocol = AbsoluteSolvationProtocol(settings=settings)

# Run the Simulation

dag = protocol.create(stateA=systemA, stateB=systemB, mapping=None)
path = pathlib.Path("./ahfe_results")
path.mkdir(exist_ok=True)
dag_results = execute_DAG(dag, scratch_basedir=path, shared_basedir=path, n_retries=3)

# Analyse the Results

protocol_results = protocol.gather([dag_results])
print(
    f"AHFE dG: {protocol_results.get_estimate()}, err {protocol_results.get_uncertainty()}"
)

outdict = {
    "estimate": protocol_results.get_estimate(),
    "uncertainty": protocol_results.get_uncertainty(),
    "protocol_result": protocol_results.to_dict(),
    "unit_results": {
        unit.key: unit.to_keyed_dict() for unit in dag_results.protocol_unit_results
    },
}

output_dir = pathlib.Path("ahfe_json")
output_dir.mkdir(exist_ok=True)
with open(output_dir / "benzene_results.json", "w") as stream:
    json.dump(outdict, stream, cls=gufe.tokenization.JSON_HANDLER.encoder)

# Absolute binding free energy of toluene to T4 lysozyme

import json
import pathlib

import gufe.tokenization
import openfe
from gufe.protocols import execute_DAG
from openfe.protocols.openmm_afe import AbsoluteBindingProtocol
from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges
from openfe.protocols.openmm_utils.omm_settings import OpenFFPartialChargeSettings
from openff.units import unit
from rdkit import Chem

# Load the Ligand

supp = Chem.SDMolSupplier("inputs/toluene.sdf", removeHs=False)
ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in supp]
print(ligands)

# Determine Charges

charge_settings = OpenFFPartialChargeSettings(
    partial_charge_method="am1bcc", off_toolkit_backend="ambertools"
)

ligands = bulk_assign_partial_charges(
    molecules=ligands,
    overwrite=False,
    method=charge_settings.partial_charge_method,
    toolkit_backend=charge_settings.off_toolkit_backend,
    generate_n_conformers=charge_settings.number_of_conformers,
    nagl_model=charge_settings.nagl_model,
    processors=1,
)

# Create Chemical Systems

solvent = openfe.SolventComponent()
protein = openfe.ProteinComponent.from_pdb_file("inputs/t4_lysozyme.pdb")  # ty: ignore[invalid-argument-type]

systemA = openfe.ChemicalSystem(  # ligand is present
    {
        "ligand": ligands[0],
        "protein": protein,
        "solvent": solvent,
    },
    name=ligands[0].name,
)

systemB = openfe.ChemicalSystem(  # ligand is fully decoupled
    {
        "protein": protein,
        "solvent": solvent,
    }
)

# Define Simulation Settings

settings = AbsoluteBindingProtocol.default_settings()
settings.protocol_repeats = 1
settings.restraint_settings.host_min_distance = 0.5 * unit.nanometer
settings.restraint_settings.host_max_distance = 1.5 * unit.nanometer
settings.engine_settings.compute_platform = "CPU"

# Create the Protocol

protocol = AbsoluteBindingProtocol(settings=settings)

# Run the Simulation

dag = protocol.create(stateA=systemA, stateB=systemB, mapping=None)

path = pathlib.Path("./abfe_results")
path.mkdir()
dag_results = execute_DAG(dag, scratch_basedir=path, shared_basedir=path, n_retries=3)

# Analyse Results

protocol_results = protocol.gather([dag_results])

print(
    f"ABFE dG: {protocol_results.get_estimate()}, err {protocol_results.get_uncertainty()}"
)

outdict = {
    "estimate": protocol_results.get_estimate(),
    "uncertainty": protocol_results.get_uncertainty(),
    "protocol_result": protocol_results.to_dict(),
    "unit_results": {
        unit.key: unit.to_keyed_dict() for unit in dag_results.protocol_unit_results
    },
}

output_dir = pathlib.Path("abfe_json")
output_dir.mkdir(exist_ok=True)
with open(output_dir / "toluene_results.json", "w") as stream:
    json.dump(outdict, stream, cls=gufe.tokenization.JSON_HANDLER.encoder)

import pathlib

import matplotlib
import openfe
from gufe.protocols.protocoldag import execute_DAG
from openfe import (
    ChemicalSystem,
    ProteinComponent,
    SmallMoleculeComponent,
    SolventComponent,
)
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from openfe.setup import LomapAtomMapper
from openfe.setup.ligand_network_planning import (
    generate_lomap_network,
    generate_minimal_spanning_network,
    generate_radial_network,
)
from openfe.utils.atommapping_network_plotting import plot_atommapping_network
from openff.units import unit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

matplotlib.use("Agg")

# Visualise Ligands

ligands_rdmol = [
    mol for mol in Chem.SDMolSupplier("inputs/tyk2_ligands.sdf", removeHs=False)
]

for ligand in ligands_rdmol:
    AllChem.Compute2DCoords(ligand)

Draw.MolsToGridImage(ligands_rdmol).save("tyk2_ligands.png")

# Load Ligands

ligands_sdf = Chem.SDMolSupplier("inputs/tyk2_ligands.sdf", removeHs=False)
ligand_mols = [SmallMoleculeComponent(sdf) for sdf in ligands_sdf]

print("Ligands:")
for ligand in ligand_mols:
    print(ligand.name)

# Ligand Atom Mapping

mapper = LomapAtomMapper()
lomap_mapping = next(
    mapper.suggest_mappings(ligand_mols[0], ligand_mols[4])
).draw_to_file("ligand_mapping.png")

# Create Ligand Network

mst_network = generate_minimal_spanning_network(
    ligands=ligand_mols,
    scorer=openfe.lomap_scorers.default_lomap_score,
    mappers=[
        LomapAtomMapper(),
    ],
)
plot_atommapping_network(mst_network).savefig("mst_network.png")

lomap_network = generate_lomap_network(
    ligands=ligand_mols,
    scorer=openfe.lomap_scorers.default_lomap_score,
    mappers=[
        LomapAtomMapper(),
    ],
)
plot_atommapping_network(lomap_network).savefig("lomap_network.png")

radial_network = generate_radial_network(
    ligands=ligand_mols[1:],
    central_ligand=ligand_mols[0],
    mappers=[
        LomapAtomMapper(),
    ],
)
plot_atommapping_network(radial_network).savefig("radial_network.png")

# Inspecting Edges

mst_edges = [edge for edge in mst_network.edges]
edge = mst_edges[1]
print("Molecule A SMILES: ", edge.componentA.smiles)
print("Molecule B SMILES: ", edge.componentB.smiles)
print("Map: ", edge.componentA_to_componentB)

with open("network_store.graphml", "w") as f:
    f.write(mst_network.to_graphml())

# Define Chemical Systems for One Edge

protein = ProteinComponent.from_pdb_file("inputs/tyk2_protein.pdb")  # ty: ignore[invalid-argument-type]
solvent = SolventComponent(
    positive_ion="Na",
    negative_ion="Cl",
    neutralize=True,
    ion_concentration=0.15 * unit.molar,
)

ejm_31_to_ejm_47 = next(
    edge for edge in mst_network.edges if edge.componentB.name == "lig_ejm_47"
)

ejm_31_complex = ChemicalSystem(
    {
        "ligand": ejm_31_to_ejm_47.componentA,
        "solvent": solvent,
        "protein": protein,
    },
    name=ejm_31_to_ejm_47.componentA.name,
)
ejm_31_solvent = ChemicalSystem(
    {
        "ligand": ejm_31_to_ejm_47.componentA,
        "solvent": solvent,
    },
    name=ejm_31_to_ejm_47.componentA.name,
)

ejm_47_complex = ChemicalSystem(
    {
        "ligand": ejm_31_to_ejm_47.componentB,
        "solvent": solvent,
        "protein": protein,
    },
    name=ejm_31_to_ejm_47.componentB.name,
)
ejm_47_solvent = ChemicalSystem(
    {
        "ligand": ejm_31_to_ejm_47.componentB,
        "solvent": solvent,
    },
    name=ejm_31_to_ejm_47.componentB.name,
)

# Define Simulation Settings

solvent_rbfe_settings = RelativeHybridTopologyProtocol.default_settings()
solvent_rbfe_settings.simulation_settings.equilibration_length = (
    10 * unit.picosecond
)  # reduce equilibration length
solvent_rbfe_settings.simulation_settings.production_length = (
    50 * unit.picosecond
)  # reduce prodution length
solvent_rbfe_settings.engine_settings.compute_platform = None

# Create the complex settings
complex_rbfe_settings = RelativeHybridTopologyProtocol.default_settings()
complex_rbfe_settings.simulation_settings.equilibration_length = (
    10 * unit.picosecond
)  # reduce equilibration length
complex_rbfe_settings.simulation_settings.production_length = (
    50 * unit.picosecond
)  # reduce equilibration length
complex_rbfe_settings.solvation_settings.solvent_padding = 1 * unit.nanometer
complex_rbfe_settings.engine_settings.compute_platform = None

print("Protocol Repeats:")
print(complex_rbfe_settings.protocol_repeats)

print("Simulation Settings:")
print(complex_rbfe_settings.simulation_settings)

print("Alchemical Settings:")
print(complex_rbfe_settings.alchemical_settings)

print("Lambda Settings:")
print(complex_rbfe_settings.lambda_settings)

print("Force Field Settings:")
print(complex_rbfe_settings.forcefield_settings)

# Create the RBFE Protocol

solvent_rbfe_protocol = RelativeHybridTopologyProtocol(settings=solvent_rbfe_settings)

complex_rbfe_protocol = RelativeHybridTopologyProtocol(settings=complex_rbfe_settings)

# Create the Transformations

transformation_complex = openfe.Transformation(
    stateA=ejm_31_complex,
    stateB=ejm_47_complex,
    mapping=ejm_31_to_ejm_47,
    protocol=complex_rbfe_protocol,
    name=f"{ejm_31_complex.name}_{ejm_47_complex.name}_complex",
)
transformation_solvent = openfe.Transformation(
    stateA=ejm_31_solvent,
    stateB=ejm_47_solvent,
    mapping=ejm_31_to_ejm_47,
    protocol=solvent_rbfe_protocol,
    name=f"{ejm_31_solvent.name}_{ejm_47_solvent.name}_solvent",
)

complex_dag = transformation_complex.create()
solvent_dag = transformation_solvent.create()

# Dump Transformations to JSON

transformation_dir = pathlib.Path("tyk2_json")
transformation_dir.mkdir(exist_ok=True)
transformation_complex.to_json(
    transformation_dir / f"{transformation_complex.name}.json"
)
transformation_solvent.to_json(
    transformation_dir / f"{transformation_solvent.name}.json"
)

# Run and Analyse

complex_path = pathlib.Path("./complex")
complex_path.mkdir()
complex_dag_results = execute_DAG(
    complex_dag, scratch_basedir=complex_path, shared_basedir=complex_path
)

solvent_path = pathlib.Path("./solvent")
solvent_path.mkdir()
solvent_dag_results = execute_DAG(
    solvent_dag, scratch_basedir=solvent_path, shared_basedir=solvent_path
)

complex_results = complex_rbfe_protocol.gather([complex_dag_results])
solvent_results = solvent_rbfe_protocol.gather([solvent_dag_results])

print(
    f"Complex dG: {complex_results.get_estimate()}, err {complex_results.get_uncertainty()}"
)
print(
    f"Solvent dG: {solvent_results.get_estimate()}, err {solvent_results.get_uncertainty()}"
)

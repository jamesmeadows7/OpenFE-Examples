import pathlib

import matplotlib
import openfe
from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges
from openfe.protocols.openmm_utils.omm_settings import OpenFFPartialChargeSettings
from openfe.utils.atommapping_network_plotting import plot_atommapping_network
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

matplotlib.use("Agg")

output_dir = pathlib.Path("ligand_network")
output_dir.mkdir(exist_ok=True)

# Visualise Ligands

ligands_rdmol = [
    mol for mol in Chem.SDMolSupplier("inputs/ligands.sdf", removeHs=False)
]
for ligand in ligands_rdmol:
    AllChem.Compute2DCoords(ligand)
Draw.MolsToGridImage(ligands_rdmol).save(output_dir / "ligands.png")

# Load Ligands

ligands_sdf = Chem.SDMolSupplier("inputs/ligands.sdf", removeHs=False)
ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in ligands_sdf]

# Assign Partial Charges

charge_settings = OpenFFPartialChargeSettings(
    partial_charge_method="am1bcc", off_toolkit_backend="ambertools"
)

charged_ligands = bulk_assign_partial_charges(
    molecules=ligands,
    overwrite=False,
    method=charge_settings.partial_charge_method,
    toolkit_backend=charge_settings.off_toolkit_backend,
    generate_n_conformers=charge_settings.number_of_conformers,
    nagl_model=charge_settings.nagl_model,
    processors=1,
)

# Create Ligand Network

mapper = openfe.LomapAtomMapper(max3d=1.0, element_change=False)
scorer = openfe.lomap_scorers.default_lomap_score
network_planner = openfe.ligand_network_planning.generate_lomap_network

ligand_network = network_planner(
    ligands=charged_ligands, mappers=[mapper], scorer=scorer
)
plot_atommapping_network(ligand_network).savefig(output_dir / "ligand_network.png")

# Inspect Edges

for edge in ligand_network.edges:
    name = f"{edge.componentA.name}_{edge.componentB.name}"
    print(f"Edge: {name}")
    print("  Molecule A SMILES: ", edge.componentA.smiles)
    print("  Molecule B SMILES: ", edge.componentB.smiles)
    print("  Score: ", edge.annotations.get("score"))
    print("  Map: ", edge.componentA_to_componentB)
    edge.draw_to_file(str(output_dir / f"mapping_{name}.png"))

# Save Network

with open(output_dir / "ligand_network.graphml", "w") as f:
    f.write(ligand_network.to_graphml())

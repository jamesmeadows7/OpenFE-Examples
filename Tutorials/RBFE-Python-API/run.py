import pathlib

import matplotlib
import openfe
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges
from openfe.protocols.openmm_utils.omm_settings import OpenFFPartialChargeSettings
from openfe.utils.atommapping_network_plotting import plot_atommapping_network
from openff.units import unit
from rdkit import Chem

matplotlib.use("Agg")

# Load Ligands

supp = Chem.SDMolSupplier("inputs/tyk2_ligands.sdf", removeHs=False)
ligands = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in supp]

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

# Create the Ligand Network

mapper = openfe.LomapAtomMapper(max3d=1.0, element_change=False)
scorer = openfe.lomap_scorers.default_lomap_score
network_planner = openfe.ligand_network_planning.generate_minimal_spanning_network

ligand_network = network_planner(
    ligands=charged_ligands, mappers=[mapper], scorer=scorer
)

plot_atommapping_network(ligand_network).savefig("ligand_network.png")

with open("ligand_network.graphml", mode="w") as f:
    f.write(ligand_network.to_graphml())

# Define Components of the Chemical Systems

solvent = openfe.SolventComponent()
protein = openfe.ProteinComponent.from_pdb_file("inputs/tyk2_protein.pdb")  # ty: ignore[invalid-argument-type]

# Create the Protocols

solvent_protocol = RelativeHybridTopologyProtocol(
    RelativeHybridTopologyProtocol.default_settings()
)

complex_settings = RelativeHybridTopologyProtocol.default_settings()
complex_settings.solvation_settings.solvent_padding = 1 * unit.nanometer
complex_protocol = RelativeHybridTopologyProtocol(complex_settings)

# Create an AlchemicalNetwork of Transformations

transformations = []
for mapping in ligand_network.edges:
    for leg in ["solvent", "complex"]:
        sysA_dict = {"ligand": mapping.componentA, "solvent": solvent}
        sysB_dict = {"ligand": mapping.componentB, "solvent": solvent}

        if leg == "complex":
            protocol = complex_protocol
            sysA_dict["protein"] = protein
            sysB_dict["protein"] = protein
        else:
            protocol = solvent_protocol

        sysA = openfe.ChemicalSystem(sysA_dict, name=f"{mapping.componentA.name}_{leg}")
        sysB = openfe.ChemicalSystem(sysB_dict, name=f"{mapping.componentB.name}_{leg}")

        prefix = "rbfe_"

        transformation = openfe.Transformation(
            stateA=sysA,
            stateB=sysB,
            mapping=mapping,
            protocol=protocol,
            name=f"{prefix}{sysA.name}_{sysB.name}",
        )
        transformations.append(transformation)

network = openfe.AlchemicalNetwork(transformations)

# Write to Disk

transformation_dir = pathlib.Path("transformations")
transformation_dir.mkdir(exist_ok=True)
for transformation in network.edges:
    transformation.to_json(transformation_dir / f"{transformation.name}.json")

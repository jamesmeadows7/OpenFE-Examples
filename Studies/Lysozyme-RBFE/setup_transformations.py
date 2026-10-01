import pathlib

import openfe
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
from openff.units import unit

# Load Ligand Network

ligand_network = openfe.LigandNetwork.from_graphml(
    pathlib.Path("ligand_network/ligand_network.graphml").read_text()
)

# Define Components of the Chemical Systems

solvent = openfe.SolventComponent()
protein = openfe.ProteinComponent.from_pdb_file(pathlib.Path("inputs/t4_lysozyme.pdb"))

# Create the Solvent Protocol

solvent_settings = RelativeHybridTopologyProtocol.default_settings()
solvent_settings.engine_settings.compute_platform = None
solvent_protocol = RelativeHybridTopologyProtocol(settings=solvent_settings)

# Create the Complex Protocol

complex_settings = RelativeHybridTopologyProtocol.default_settings()
complex_settings.solvation_settings.solvent_padding = 1 * unit.nanometer
complex_settings.engine_settings.compute_platform = None
complex_protocol = RelativeHybridTopologyProtocol(settings=complex_settings)

# Create an AlchemicalNetwork of Transformations

prefix = "rbfe_"
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

        transformation = openfe.Transformation(
            stateA=sysA,
            stateB=sysB,
            mapping=mapping,
            protocol=protocol,
            name=f"{prefix}{sysA.name}_{sysB.name}",
        )
        transformations.append(transformation)

network = openfe.AlchemicalNetwork(transformations)

# Write Transformations to Disk

transformation_dir = pathlib.Path("transformations")
transformation_dir.mkdir(exist_ok=True)

for transformation in network.edges:
    transformation.to_json(transformation_dir / f"{transformation.name}.json")

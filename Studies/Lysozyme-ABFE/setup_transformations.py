import pathlib

import openfe
from openfe.protocols.openmm_afe import AbsoluteBindingProtocol
from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges
from openfe.protocols.openmm_utils.omm_settings import OpenFFPartialChargeSettings
from openff.units import unit
from rdkit import Chem

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

# Define Components of the Chemical Systems

solvent = openfe.SolventComponent()
protein = openfe.ProteinComponent.from_pdb_file(pathlib.Path("inputs/t4_lysozyme.pdb"))


# Create the Protocol

settings = AbsoluteBindingProtocol.default_settings()
settings.protocol_repeats = 1
settings.restraint_settings.host_min_distance = 0.5 * unit.nanometer
settings.restraint_settings.host_max_distance = 1.5 * unit.nanometer
settings.engine_settings.compute_platform = None
protocol = AbsoluteBindingProtocol(settings=settings)

# Create the Transformations

prefix = "abfe_"
transformations = []
for ligand in charged_ligands:
    sysA = openfe.ChemicalSystem(  # ligand present
        {
            "ligand": ligand,
            "protein": protein,
            "solvent": solvent,
        },
        name=ligand.name,
    )

    sysB = openfe.ChemicalSystem(  # ligand fully decoupled
        {
            "protein": protein,
            "solvent": solvent,
        }
    )

    transformation = openfe.Transformation(
        stateA=sysA,
        stateB=sysB,
        mapping=None,
        protocol=protocol,
        name=f"{prefix}{ligand.name}",
    )
    transformations.append(transformation)

# Write Transformations to Disk

transformation_dir = pathlib.Path("transformations")
transformation_dir.mkdir(exist_ok=True)

for transformation in transformations:
    transformation.to_json(transformation_dir / f"{transformation.name}.json")

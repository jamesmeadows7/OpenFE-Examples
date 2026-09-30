# OpenFE Examples

Learning [OpenFE](https://docs.openfree.energy/en/stable/index.html) through a series of tutorials.

Techniques covered:

- Atom Mapping (LOMAP)
- Ligand Network Planning (Minimal Spanning Tree, LOMAP, Radial)
- Deriving AM1-BCC Partial Charges
- Relative Binding Free Energy (RBFE) Calculations
- Absolute Binding Free Energy (ABFE) Calculations
- Absolute Hydration Free Energy Calculations

| # | Tutorial | Description | Script |
|---|----------|-------------|--------|
| 1 | [OpenFE Showcase](https://docs.openfree.energy/en/stable/tutorials/showcase_notebook.html) | RBFE workflow from atom mapping and ligand networks to running transformation calculations | [run.py](Tutorials/OpenFE-Showcase/run.py) |
| 2 | [MD Protocol](https://docs.openfree.energy/en/stable/tutorials/md_tutorial.html) | Plain MD simulation of benzene bound to T4 lysozyme | [run.py](Tutorials/MD-Protocol/run.py) |
| 3 | [RBFE with the Python API](https://docs.openfree.energy/en/stable/tutorials/rbfe_python_tutorial.html) | RBFE campaign for a set of ligands bound to TYK2 protein | [run.py](Tutorials/RBFE-Python-API/run.py) |
| 4 | [ABFE Protocol](https://docs.openfree.energy/en/stable/tutorials/abfe_tutorial.html) | Absolute binding free energy of toluene to T4 lysozyme | [run.py](Tutorials/ABFE-Protocol/run.py) |
| 5 | [Solvation Free Energy Protocol](https://docs.openfree.energy/en/stable/tutorials/ahfe_tutorial.html) | Absolute hydration free energy of benzene in water | [run.py](Tutorials/Solvation-Free-Energy-Protocol/run.py) |

## Setup

```bash
conda env create -f environment.yml
conda activate openfe
```

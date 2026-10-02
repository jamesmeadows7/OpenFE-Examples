#!/bin/bash
#SBATCH --job-name=t4l_abfe
#SBATCH --partition=cuda
#SBATCH --time=2-00:00:00
#SBATCH --gres=gpu:h200_nvl_1g.18gb:1
#SBATCH --array=0-2

source /nobackup/drrj44/miniforge3/bin/activate openfe

transformations=(transformations/*.json)
transformation=${transformations[$SLURM_ARRAY_TASK_ID]}
name=$(basename "$transformation" .json)

echo "Running $name on $(hostname)"
nvidia-smi
python -m openmm.testInstallation

mkdir -p results simulations/"$name"
openfe quickrun "$transformation" -o results/"$name".json -d simulations/"$name"

#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=28G

srun "$@"

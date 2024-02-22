#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=nf_default
#SBATCH --error=slurm-%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kk4764@nyu.edu

module purge
module load nextflow/23.04.1

nextflow run nf-core/sarek -r 3.4.0 --input ./samplesheet.csv --outdir res2 -profile nyu_hpc

echo _ESTATUS_ [ nextflow run ]: $?
echo _END_ [ test_sarek.slurm ]: $(date)


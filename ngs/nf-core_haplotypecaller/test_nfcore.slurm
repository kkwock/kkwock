#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=nf_template
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kk4764@nyu.edu

module purge
module load nextflow/23.04.1

sarek_folder=/scratch/work/courses/BI7653/hw4.2024/sarek
kk_folder=/scratch/kk4764/nextflow/ngs

nextflow run nf-core/sarek -r 3.4.0 --input ./samplesheet.csv --outdir res3_10M -profile nyu_hpc -c ${sarek_folder}/sarek.nyu_ngs.config

echo _ESTATUS_ [ nextflow run ]: $?
echo _END_ [ test_sarek.slurm ]: $(date)


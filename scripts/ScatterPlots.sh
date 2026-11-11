#!/bin/bash

#PBS -S /bin/sh
#PBS -N BCC_scatter
#PBS -l select=1:ncpus=8:mem=500G
#PBS -l walltime=120:00:00
#PBS -q longq

source /home/vfama/miniconda3/bin/activate /home/vfama/miniconda3/envs/genomicfeatures

cd /projects/CGS_shared/vfama/BRIGHT_PROJECT/scripts

Rscript Stoichiometry_reduction_scatterplots.R

conda deactivate
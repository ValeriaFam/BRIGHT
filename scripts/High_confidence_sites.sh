#!/bin/bash

#PBS -S /bin/sh
#PBS -N High_conf_sites
#PBS -l select=1:ncpus=4:mem=250G
#PBS -l walltime=12:00:00
#PBS -q longq

source /home/vfama/miniconda3/bin/activate /home/vfama/miniconda3/envs/genomicfeatures2

cd /projects/CGS_shared/vfama/BRIGHT_PROJECT/scripts

Rscript High_confidence_sites_and_DEGs.R

conda deactivate
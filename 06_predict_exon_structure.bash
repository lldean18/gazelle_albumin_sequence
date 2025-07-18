#!/bin/bash
# Laura Dean
# 16/7/25
# script written for running on the UoN HPC Ada

# script to predict the exon strucuture of the albumin protein based on albumin sequences from
# closely related / other species and the DNA sequence of the region containing the gene in
# Gazella cuvieri

#SBATCH --job-name=predict_exon_structure
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/gazelle
albumin_seqs=( /gpfs01/home/mbzlld/data/gazelle/bams/clean_bams/*_albumin_region.fasta )

cd $wkdir

# activate software
source $HOME/.bash_profile
#conda create --name exonerate bioconda::exonerate -y
conda activate exonerate

for seq in ${albumin_seqs[@]}; do
# predict the exon structure
exonerate --model protein2genome \
          --showtargetgff yes \
          --showalignment no \
          --query albumin_protein.fasta \
          --target $seq > ${seq%.*}_exonerate.out
done

# deactivate software
conda deactivate


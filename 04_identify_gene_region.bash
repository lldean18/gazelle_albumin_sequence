#!/bin/bash
# Laura Dean
# 14/7/25
# script written for running on the UoN HPC Ada

#SBATCH --job-name=identify_gene_region
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out



wkdir=/gpfs01/home/mbzlld/data/gazelle
cd $wkdir


# activate software
source $HOME/.bash_profile
#conda create --name blast bioconda::blast -y
conda activate blast



# blast albumin protein sequences against consensus sequences
tblastn -query $wkdir/albumin_protein.fasta \
        -subject $wkdir/bams/clean_bams/ERR7570036_consensus.fasta \
        -out $wkdir/albumin_tblastn.out \
        -outfmt 6 \
        -evalue 1e-5 \
        -num_threads 4



# deactivate software
conda deactivate


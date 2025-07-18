#!/bin/bash
# Laura Dean
# 16/7/25
# script written for running on the UoN HPC Ada

# script to extract the region containing the albumin gene from Gazella cuvieri
# consensus sequences for each individual

#SBATCH --job-name=extract_gene_region
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10g
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/gazelle
ref=$wkdir/CAKJTW01.fasta
consensus_seqs=( /gpfs01/home/mbzlld/data/gazelle/bams/clean_bams/*_consensus.fasta )

cd $wkdir

# activate software
source $HOME/.bash_profile
#conda create --name samtools1.22 bioconda::samtools=1.22 bioconda::bcftools=1.22 -y
#conda install seqtk -y # code doesn't actually need this any more
conda activate samtools1.22

for seq in ${consensus_seqs[@]}; do
# extract the gene region
samtools faidx --threads 15 $seq
samtools faidx --threads 15 $seq 'ENA|CAKJTW010000001|CAKJTW010000001.1:16846815-16909313' > ${seq%.*}_albumin_region.fasta
done

# deactivate software
conda deactivate


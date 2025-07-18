#!/bin/bash
# Laura Dean
# 14/7/25
# script written for running on the UoN HPC Ada

# script to identify the region of the Gazella cuvieri genome that contains the
# albumin gene based on a blast search of closely related protein sequences.

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
        -evalue 1e-10 \
        -num_threads 4


# deactivate software
conda deactivate


# The second column of the ouput gives the chr name
# The 9th and 10th columns give the start and end of the match
# The sheep albumin protein is the best match (unsurprising since it is the most closely related)
# it must be on the reverse strand bc the start and end are the wrong way around
# full stretch of sheep match is 16856815 to 16899313 of ENA|CAKJTW010000001|CAKJTW010000001.1 which spans 42,498 bases
# will take 10,000 bases extra from each side so 16846815 to 16909313








#!/bin/bash
# Laura Dean
# 2/7/25
# script written for running on the UoN HPC Ada

# script to generate a concensus sequence from reads aligned to a reference

#SBATCH --job-name=generate_consensus
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10g
#SBATCH --time=50:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/gazelle
ref=$wkdir/CAKJTW01.fasta
bams=( /gpfs01/home/mbzlld/data/gazelle/bams/clean_bams/*.bam )

cd $wkdir

# activate software
source $HOME/.bash_profile
#conda create --name samtools1.22 bioconda::samtools=1.22 bioconda::bcftools=1.22 -y
#conda install seqtk -y
conda activate samtools1.22

for bam in ${bams[@]}; do
# generate the concensus sequence
bcftools mpileup --fasta-ref $ref $bam | bcftools call --consensus-caller -Oz -o calls.vcf.gz
bcftools index calls.vcf.gz
bcftools consensus --fasta-ref $ref calls.vcf.gz > ${bam%.*}_consensus.fasta
rm calls.vcf.gz calls.vcf.gz.csi
done

#bcftools mpileup -Ou --fasta-ref $ref $bam | bcftools call --consensus-caller - | vcfutils.pl vcf2fq - > ${bam%.*}_consensus.fq
#seqtk seq -a ${bam%.*}_consensus.fq > ${bam%.*}_consensus.fa
##rm ${bam%.*}_consensus.fq


# deactivate software
conda deactivate


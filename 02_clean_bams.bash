#!/bin/bash
# Laura Dean
# 11/7/25
# script written for running on the UoN HPC Ada

# script to clean and filter bam files
# files contain Gazella cuvieri reads algned to the Nangar dama reference

#SBATCH --job-name=clean_bams
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=19
#SBATCH --mem=30g
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out
#SBATCH --array=1-10


# set variables
wkdir=/gpfs01/home/mbzlld/data/gazelle/bams
config=/gpfs01/home/mbzlld/code_and_scripts/config_files/gazelle_array_config.txt
individual=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)


# load software
module load samtools-uoneasy/1.18-GCC-12.3.0
module load bcftools-uoneasy/1.18-GCC-13.2.0



# print info to slurm output
echo "This is array task ${SLURM_ARRAY_TASK_ID}, cleaning individual $individual,
cleaned output BAM files will be written to the folder $wkdir/clean_bams"



# clean bams
# Remove unmapped reads and do quality filtering
# -q mapping quality greater than or equal to 40
# -f include reads mapped in a propper pair
# -F Only include reads which are not read unmapped or mate unmapped
samtools view \
--threads 19 \
-q 40 \
-f 2 \
-F 4 \
-b $wkdir/raw_bams/${individual}_raw.bam |
# Mark and remove duplicate reads (the -r flag removes the duplicate reads, removing it will mark dups only)
samtools markdup -r --threads 19 - $wkdir/clean_bams/$individual.bam

# index the final BAM files
samtools index -@ 19 $wkdir/clean_bams/$individual.bam

# check the mapping
echo "after cleaning and filtering the final mapping success was:"
samtools flagstat $wkdir/clean_bams/$individual.bam




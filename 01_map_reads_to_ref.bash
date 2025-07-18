#!/bin/bash
# Laura Dean
# 10/7/25
# script written for running on the UoN HPC Ada

# script to align raw sequencing reads for Gazella cuvieri (downloaded from ENA)
# to the Nangar dama reference genome (downloaded from ENA)
# this is the most closely related available reference genome

#SBATCH --job-name=map_gazelle_reads
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=19
#SBATCH --mem=30g
#SBATCH --time=35:00:00
#SBATCH --array=1-10
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out





########################################################
# SET UP YOUR ENVIRONMENT AND SPECIFY DATA INFORMATION #
########################################################

# specify your config file
config=/gpfs01/home/mbzlld/code_and_scripts/config_files/gazelle_array_config.txt

# extract the sample name for the current slurm task ID
individual=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# set the input data location
in_filepath=/gpfs01/home/mbzlld/data/gazelle

# set the output filepath where you want to put the bams
out_filepath=/gpfs01/home/mbzlld/data/gazelle

#specify the reference genome you will use
reference_genome=/gpfs01/home/mbzlld/data/gazelle/CAKJTW01.fasta

echo "This is array task ${SLURM_ARRAY_TASK_ID}, aligning reads for individual $individual, from the input filepath $in_filepath
output BAM files will be written to the folder $out_filepath/bams"

# load the necessary modules
module load bwa-uoneasy/0.7.17-GCCcore-12.3.0
module load samtools-uoneasy/1.18-GCC-12.3.0
module load bcftools-uoneasy/1.18-GCC-13.2.0

# Index the reference genome you want to use (only need to do this once)
#bwa index $reference_genome

# name the paired files for each individual
R1filename=$""$in_filepath"/"$individual"_1.fastq.gz"
R2filename=$""$in_filepath"/"$individual"_2.fastq.gz"
echo ""$R1filename" 
"$R2filename" 
are the filenames being worked on"

# Extract the header line of the fastq file
file_info=$(zcat $R1filename | head -n 1)

# Save the pieces of information you need as variables
flowcell_ID=$(cut -d ":" -f3 <<< "$file_info")
lane_no=$(cut -d ":" -f4 <<< "$file_info")
sample_barcode=$(cut -d ":" -f10 <<< "$file_info")

# store the read group information for the ID and PU fields as variables from the individual ones you just created
PU=$flowcell_ID.$lane_no.$sample_barcode
ID=$flowcell_ID.$lane_no

#######################################
# ALIGN READS TO THE REFERENCE GENOME #
#######################################

# make a directory in which to put the BAM files if it doesn't already exist
mkdir -p $out_filepath/bams
# within this directory make two directories to contain the cleaned and raw bam files
mkdir -p $out_filepath/bams/raw_bams

###### Align the reads to the reference genome using bwa mem ######
# BWA MEM command explanation:
# -t = number of threads
# -M = Mark shorter split hits as secondary (for Picard compatibility)
# -R = Read group header line - gives identifiers for indivs and batches (lanes). Allows merging files which contain data for same indivs. Need unique ID and SM
#      ID: the read group identifier - the flowcell name followed by lane number this info is in the first line of the fastq file 
#      PU: the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
#      SM: the name of the individual
#      PL: the platform used for sequencing (ILLUMINA)
#      LB: Unique ID for the library prep (this is the same as the individual ID unless you pooled DNA before lib prep)
# -o /path/to/output/sam/file
# /path/to/reference/genome
# /path/to/forward/reads
# /path/to/reverse/reads
bwa mem \
-t 19 \
-M \
-R "@RG\tID:"$ID"\tSM:"$individual"\tPL:ILLUMINA\tLB:"$individual"\tPU:"$PU"" \
$reference_genome \
$R1filename \
$R2filename |
# mark the duplicate reads then
# Sort the SAM files
samtools fixmate --threads 19 -m -O BAM - - |
samtools sort --threads 19 -o $out_filepath/bams/raw_bams/$individual\_raw.bam

# Index the BAM files
samtools index -@ 19 $out_filepath/bams/raw_bams/$individual\_raw.bam

# Generate info  - look at how well the reads mapped
echo "the raw reads mapped with the following success:"
samtools flagstat --threads 19 $out_filepath/bams/raw_bams/$individual\_raw.bam


# to check that the file header has been printed correctly
#module load samtools-uoneasy/1.12-GCC-9.3.0
#samtools view -H /gpfs01/home/mbzlld/data/Ponds/Ponds2022/bams/Ponds_22A009.bam | grep '^@RG'




#!/bin/bash
# Laura Dean
# 17/7/25
# script written for running on the UoN HPC Ada

# script to process a predicted annotation of a gene to protein sequence

#SBATCH --job-name=convert_annotation_to_protein
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/gazelle
predicted_annotations=( /gpfs01/home/mbzlld/data/gazelle/bams/clean_bams/*_exonerate.out )


cd $wkdir


# activate software
source $HOME/.bash_profile
conda activate bedtools 

for annotation in ${predicted_annotations[@]}; do
# === extract sheep gff section ===
sed -n '/^vulgar: tr|W5PWE9|W5PWE9_SHEEP/,/^# --- END OF GFF DUMP ---/p' $annotation > ${annotation%.*}_sheep.gff
# === extract exons ===
grep -P "\texon\t" ${annotation%.*}_sheep.gff > ${annotation%.*}_sheep_exons.gff
# === get the strand ===
STRAND=$(awk '{print $7}' ${annotation%.*}_sheep_exons.gff | uniq)
# === convert exons to fasta ===
bedtools getfasta -fi ${annotation%_*}.fasta -bed ${annotation%.*}_sheep_exons.gff -s -name > ${annotation%.*}_sheep_exons.fasta
# === merge exons into single CDS FASTA ===
# Remove headers, flatten, add new header
seqkit seq ${annotation%.*}_sheep_exons.fasta | grep -v ">" | tr -d '\n' > cds.tmp
echo ">Albumin_cds" > ${annotation%_*}_cds.fasta
cat cds.tmp >> ${annotation%_*}_cds.fasta
rm cds.tmp
# === reverse-complement if minus strand ===
#if [ "$STRAND" = "-" ]; then
#  echo "Reverse-complementing because gene is on minus strand..."
#  seqkit seq -r -p -t DNA ${annotation%_*}_cds.fasta > ${annotation%_*}_cds_rc.fasta
#  mv ${annotation%_*}_cds_rc.fasta ${annotation%_*}_cds.fasta
#fi
# === translate cds to protein ===
seqkit translate ${annotation%_*}_cds.fasta > ${annotation%_*}_protein.fasta

#rm ${annotation%.*}_sheep.gff
#rm ${annotation%.*}_sheep_exons.gff
#rm ${annotation%.*}_sheep_exons.bed
#rm ${annotation%.*}_sheep_exons_sorted.bed
#rm ${annotation%.*}_sheep_exons.fasta
done

# deactivate software
conda deactivate


#!/bin/bash
# Laura Dean
# 18/7/25
# script written for running on the UoN HPC Ada

# script to align the individual gazelle protein sequences and generate a single consensus
# sequence for the species

#SBATCH --job-name=align_proteins
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=30g
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/mbzlld/code_and_scripts/slurm_out_scripts/slurm-%x-%j.out

# set variables
wkdir=/gpfs01/home/mbzlld/data/gazelle
cd $wkdir


# combine fasta files
cat $wkdir/bams/clean_bams/*_protein.fasta > $wkdir/Gazelle_albumin_protein_seqs.fasta

# align the sequences by pasting the fasta into https://www.genome.jp/tools-bin/clustalw
# on inspecting the alignment, sequences are the same apart from the occasional X and a big chunk 
# that doesn't match up in ERR7570036 which is probably a result of poor sequencing in that region 
# in that one sample.

# then generate a consensus protein sequence for the Gazelle by pasting the alignment into 
# https://www.ebi.ac.uk/jdispatcher/msa/emboss_cons?stype=protein

# take the consensus sequence generated above and put into a file
echo ">Gazella_cuvieri_albumin
MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEENFQGLVLIAFSQYLQQCPF
DDHVKLVKEVTEFAKTCVADESHAGCDKSLHTIFGDELCKVATLRETYGDMADCCEKQEP
ERNECFLKHKDDSPDLPKLKPEPDTLCAEFKADEKKFWGKYLYEVARRHPYFYAPELLYY
ANKYNGVFQECCEAEDKGACLLPKIETMREKVLASSARQRLKCASIQKFGERALKAWSVA
RLSQKFPKADFTEVTKLVEDVTKVHKECCHGDLLECADDRADLAKYMCDHQDAISGKLKE
CCDKPLLEKSHCLAEVDEDTMPENLAPLTAAFAEDKEVCKNYQEAKDIFLGSFLYEYSRR
HPDYAVSVLLRLAKEYEATLEDCCAKDDPHACYATVFDKLKHLVDEPQDLIKKNCELFEK
HGEYGFQNALIVRYTRKVPQVSTPTLVEISRSLGKVGTKCCTQPESKRMPCTEDYLSLIL
NRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTLDETYVPKPFDDKLFTFHADICTLP
DTEKQIKKQTCSALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEGCFALEGPKL
VASTQAALA" > $wkdir/Gazella_cuvieri_albumin.fasta






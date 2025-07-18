# gazelle_albumin_sequence

Repo for scripts used to determine the protein sequence of Albumin in Gazella cuvieri.

This genus has no published assemblies or annotations! (whaat!?) so the pipeline involves:

1. Mapping illumina short-read sequencing data to the most closely related available reference genome (Nangar dama).

2. Cleaning up and filtering those mapped reads.

3. Generating consensus assemblies for each Gazella cuvieri individual based on the mapped reads. This generates pseudo reference assemblies for each individual.

4. Identifying the region of the Gazella cuvieri pseudo-genomes that contains the Albumin coding region by blasting protein sequences from other species against the pseudo-genomes.

5. Extracting the gene region (and some buffer surrounding it) that codes for the Albumin protein from each pseudo-genome.

6. Using the coding region and known protein sequences to predict the exon strucutre of the Albumin protein in Gazella cuvieri.

7. Extracting and concatenating those exon sequences and converting them to protein sequences.

8. Aligning the predicted Albumin protein sequences for each individual and extracting the consensus sequence for the species from the alignment. All individuals were largely the same at this point apart from a few ambiguities in single individuals here and there, and one individual with a chunk that was totally off towards the end. This suggests the process worked well and the final consensus sequence is likey reasonably accurate.




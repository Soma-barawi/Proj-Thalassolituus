#!/bin/bash
#SBATCH --account=rrg-rbeiko
#SBATCH --job-name=nif-trees
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --mail-user=sm871954@dal.ca
#SBATCH --mail-type=END,FAIL

# Phylogenetic analysis of the 6 nif genes and their corresponding Nif proteins
# Diaiden was used to annotate and extract nifH, nifD, nifK, nifE, nifN and nifB

mkdir -p nif-alignments nif-trees concatenated-alignments

# Align the nucleotide (.fna) and amino acid (.faa) sequences with MAFFT
for gene in nifH nifD nifK nifE nifN nifB; do
    mafft --auto --thread 20 "${gene}.fna" > "nif-alignments/${gene}_aligned.fna"
    mafft --auto --thread 20 "${gene}.faa" > "nif-alignments/${gene}_aligned.faa"
done

# Check all protein alignments
# Run IQ-TREE
for gene in nifH nifD nifK nifE nifN nifB; do
    mkdir -p "nif-trees/${gene}-AA-NT-trees/AA-tree"
    mkdir -p "nif-trees/${gene}-AA-NT-trees/NT-tree"

    iqtree3 -s "nif-alignments/${gene}_aligned.faa" -m MFP -B 1000 --alrt 1000 -T 20 --prefix "nif-trees/${gene}-AA-NT-trees/AA-tree/${gene}"
    iqtree3 -s "nif-alignments/${gene}_aligned.fna" -m MFP -B 1000 --alrt 1000 -T 20 --prefix "nif-trees/${gene}-AA-NT-trees/NT-tree/${gene}"
done

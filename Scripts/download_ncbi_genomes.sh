#!/bin/bash

# NCBI Datasets
# Use datasets to download biological sequence data across all domains of life from NCBI
# Use dataformat to convert metadata from JSON Lines format to other formats

# Initial Oceanospirillaceae genome dataset 
datasets download genome taxon "Oceanospirillaceae" --include genome,protein --filename NCBI-Oceanospirillaceae.zip
dataformat excel genome --inputfile NCBI-Oceanospirillaceae/ncbi_dataset/data/assembly_data_report.jsonl --outputfile NCBI-Oceanospirillaceae-metadata.xlsx

mkdir -p NCBI-Oceanospirillaceae
unzip -q NCBI-Oceanospirillaceae.zip -d NCBI-Oceanospirillaceae

# Another useful python script for downloading genomes from RefSeq/ Genbank
ncbi-genome-download --genera genomes.txt --formats fasta --assembly-levels complete,chromosome,scaffold,contig bacteria
ncbi-genome-download bacteria --assembly-accessions accessions.txt --output-folder genomes_downloaded 

#!/bin/bash

# Run ARETE on all Thalassolituus genomes

# Because I've separated RefSeq and Genbank genomes
refseq_dir="/lustre06/project/6070434/somab/Projects/ARETE/Thalassolituus/Oceanospirillaceae/all_refseq/all_fna_gz"
genbank_dir="/lustre06/project/6070434/somab/Projects/ARETE/Thalassolituus/Oceanospirillaceae/all_genbank/all_fna_gz"

# Removing .1 from the filenames because ARETE is complaining
for folder in "$refseq_dir" "$genbank_dir"; do
    for file in "$folder"/*.1.fna; do
        [ -e "$file" ] || continue
        mv -n "$file" "${file%.1.fna}.fna"
    done
done

# Make the ARETE sample sheet
# GCF accessions are from RefSeq; GCA_ are from GenBank
awk -v refseq="$refseq_dir" -v genbank="$genbank_dir" '{
    accession=$1
    sub(/\.[0-9]+$/, "", accession)

    if (accession ~ /^GCF_/) {
        print accession "," refseq "/" accession ".fna"
    } else {
        print accession "," genbank "/" accession ".fna"
    }
}' accessions.txt > samplesheet.csv

# Run ARETE
nextflow run Application-ARETE/ \
    --input_sample_table samplesheet.csv \
    --db_cache /home/somab/projects/rrg-rbeiko/somab/Projects/ARETE/Application-ARETE/dbcache \
    --bakta_db /home/somab/projects/rrg-rbeiko/somab/Projects/ARETE/Application-ARETE/dbcache/baktadb/db-light/ \
    --annotation_tools 'mobsuite,rgi,cazy,vfdb,iceberg,bacmet,islandpath,phispy,integronfinder,report' \
    --outdir all_thalasso_output \
    --run_recombination \
    --skip_poppunk \
    --run_evolccm \
    --run_rspr \
    -entry annotation \
    -profile singularity,narval \
    -resume

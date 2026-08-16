#!/bin/bash
# Screening the 25 putative diazotrophic genomes for mobile genetic elements

mkdir -p hmmsearch_results geNomad_output MGE_output hmm_profiles HMM_hits

# Search each proteome against the proMGE recombinase HMM database
hmmpress -f proMGE_all.hmm

for protein in Proteins/*.faa; do
    name=$(basename "$protein" .faa)
    hmmscan --cpu 8 --domtblout "hmmsearch_results/${name}.domtbl" --tblout "hmmsearch_results/${name}.tbl" proMGE_all.hmm "$protein" > "hmmsearch_results/${name}.txt"
done

# Run geNomad to identify plasmid and virus-associated sequences
for genome in diazotroph_genomes/*.fna; do
    name=$(basename "$genome" .fna)
    mkdir -p "geNomad_output/$name"
    apptainer run --bind "$(pwd):/app" genomad.sif end-to-end --cleanup --splits 8 "/app/$genome" "/app/geNomad_output/$name" /app/genomad_db
done

# Run MobileElementFinder 
for genome in diazotroph_genomes/*.fna; do
    name=$(basename "$genome" .fna)
    mefinder find --contig "$genome" --threads 4 "MGE_output/${name}_mge"
done

# Update and build the individual HMM profiles from the aligned proMGE protein sequences
for alignment in *.aln; do
    name=$(basename "$alignment" .aln)
    hmmbuild "hmm_profiles/${name}.hmm" "$alignment"
done

cat hmm_profiles/*.hmm > all_MGEs.hmm
hmmpress -f all_MGEs.hmm

# Search each proteome against the updated HMM database
for protein in Proteins/*.faa; do
    name=$(basename "$protein" .faa)
    hmmscan --cpu 4 -E 1e-5 --domtblout "HMM_hits/${name}.domtbl" --tblout "HMM_hits/${name}.tbl" all_MGEs.hmm "$protein" > "HMM_hits/${name}.txt"
done

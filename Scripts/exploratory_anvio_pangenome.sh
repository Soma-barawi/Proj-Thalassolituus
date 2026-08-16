#!/bin/bash

# Only for exploratory purposes; building a pangenome and phylogenomic tree in Anvi'O using the 25 putative diazotrophic genomes 
# Run locally in a conda environment

mkdir -p genomes contigs-db Thalasso

# Simplify the contig names so they work with anvi'o
for fasta in raw-genomes/*.fna; do
    name=$(basename "$fasta" .fna)
    anvi-script-reformat-fasta "$fasta" -o "genomes/${name}.fna" --simplify-names
done

# Make an anvi'o contigs database for each genome
for fasta in genomes/*.fna; do
    name=$(basename "$fasta" .fna)
    anvi-gen-contigs-database -f "$fasta" --project-name "$name" -o "contigs-db/${name}.db" -T 8 --force-overwrite
done

# Add all the functional, taxonomic and tRNA annotations to each contigs database
for db in contigs-db/*.db; do
    anvi-run-ncbi-cogs -c "$db" -T 8
    anvi-run-hmms -c "$db" -T 8
    anvi-run-kegg-kofams -c "$db" -T 8
    anvi-run-pfams -c "$db" -T 8
    anvi-run-scg-taxonomy -c "$db" -T 8
    anvi-run-trna-taxonomy -c "$db" -T 8
done

# external_genomes.txt contains the genome names and paths to their contigs databases
anvi-gen-genomes-storage -e external_genomes.txt -o GENOMES.db

# Build the pangenome
anvi-pan-genome -g GENOMES.db --project-name Thalassolituus_Pan --output-dir Thalasso --num-threads 8 --minbit 0.5 --mcl-inflation 10 --use-ncbi-blast
pan_db="Thalasso/Thalassolituus_Pan-PAN.db"

# Extract and concatenate the single copy gene clusters found once in each of the 25 genomes, then align with MAFFT
anvi-get-sequences-for-gene-clusters -p "$pan_db" -g GENOMES.db --min-num-genomes-gene-cluster-occurs 25 --max-num-genes-from-each-genome 1 --concatenate-gene-clusters --output-file Thalassolituus-SCGs.fa

# Align, trim, build the tree 
mafft --auto --reorder Thalassolituus-SCGs.fa > Thalassolituus-SCGs-aligned.fasta
trimal -in Thalassolituus-SCGs-aligned.fasta -out Thalassolituus-SCGs-trimmed.fa -gt 0.50
iqtree -s Thalassolituus-SCGs-trimmed.fa -nt 8 -m LG -bb 1000

# Make the layer-order file used to import the tree into the anvi'o pangenome
echo -e "item_name\tdata_type\tdata_value" > Thalassolituus-phylogenomic-layer-order.txt
echo -e "SCGs_Tree\tnewick\t$(cat Thalassolituus-SCGs-trimmed.fa.contree)" >> Thalassolituus-phylogenomic-layer-order.txt

# Import the phylogenomic tree and view in the interactive interface
anvi-import-misc-data -p "$pan_db" -t layer_orders Thalassolituus-phylogenomic-layer-order.txt
anvi-display-pan -g GENOMES.db -p "$pan_db"

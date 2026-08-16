# NCBI Datasets
# Use datasets to download biological sequence data across all domains of life from NCBI
# Use dataformat to convert metadata from JSON Lines format to other formats

# Initial Oceanospirillaceae genome dataset 
datasets download genome taxon "Oceanospirillaceae" --include genome,protein --filename NCBI-Oceanospirillaceae.zip
dataformat excel genome --inputfile NCBI-Oceanospirillaceae/ncbi_dataset/data/assembly_data_report.jsonl --outputfile NCBI-Oceanospirillaceae-metadata.xlsx

mkdir -p NCBI-Oceanospirillaceae
unzip -q NCBI-Oceanospirillaceae.zip -d NCBI-Oceanospirillaceae

# Another useful python script for downloading genomes from RefSeq/ Genbank

conda create -n ncbi_env python=3.9
conda update -n base -c conda-forge conda
conda activate ncbi_env
pip install ncbi-genome-download

ncbi-genome-download --genera Thalassolituus --section refseq --formats fasta, protein-fasta bacteria
ncbi-genome-download --genera genomes.txt --formats fasta --assembly-levels complete,chromosome,scaffold,contig bacteria
ncbi-genome-download bacteria --assembly-accessions accessions.txt --output-folder genomes_downloaded 

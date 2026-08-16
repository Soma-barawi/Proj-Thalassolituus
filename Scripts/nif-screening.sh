# Run locally
# Diaiden uses Prodigal to predict proteins and DIAMOND to identify nitrogen-fixation genes
# Screening the 74 Thalassolituus genomes for all 6 nif genes

# One-time installation
# conda create -n diaiden -y
# conda activate diaiden
# mamba install r-base=4.4.0 r-dplyr=1.1.4 prodigal=2.6.3 diamond=2.1.6 -y
# git clone https://github.com/jchenek/Diaiden.git

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate diaiden

genomes="/Users/somabarawi/Bioinformatics/diaiden/Diaiden/Thalassolituus_genomes"
diaiden="/Users/somabarawi/Bioinformatics/diaiden/Diaiden"

perl "$diaiden/Diaiden.pl" -i "$genomes" -p "$diaiden" -c 3 -b 3
# nif_anno_full_info.tsv contains all nif annotations
# nifh/nifd/nifk/nife/nifn/nifb .fna and .faa files contain the extracted nucleotide and protein sequences

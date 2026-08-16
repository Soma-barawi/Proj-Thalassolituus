#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --job-name=blastp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --account=rrg-rbeiko
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=sm871954@dal.ca

# Search the Thalassolituus proteins against Swiss-Prot and NCBI NR
module load StdEnv/2023 gcc/12.3
module load blast+/2.14.1

proteins="/home/somab/projects/rrg-rbeiko/somab/Projects/ARETE/Thalassolituus/thalassolituus_proteins"
db="/home/somab/projects/rrg-rbeiko/somab/Projects/ARETE/Thalassolituus/blast/swissprot_db"
output="/home/somab/projects/rrg-rbeiko/somab/Projects/ARETE/Thalassolituus/blast/swissprot_results"
mkdir -p "$output"

for protein_file in "$proteins"/*.faa; do
    name=$(basename "$protein_file" .faa)
    blastp -query "$protein_file" -db "$db" -out "$output/${name}_top_hits.txt" -max_target_seqs 1 -evalue 1e-25 -max_hsps 1 -outfmt "6 qseqid bitscore evalue length qlen qcovs pident sseqid sgi sacc staxids sscinames scomnames stitle sseq" -num_threads "$SLURM_CPUS_PER_TASK"
done

# NCBI non-redundant protein database
nr_queries="/home/somab/links/projects/rrg-rbeiko/somab/Proj-Thalassolituus/All_Thalassolituus_analysis/BLAST_analyses/nif_queries"
nr_db="/cvmfs/bio.data.computecanada.ca/content/databases/Core/blast_dbs/2025_06_21/nr/nr"
nr_output="/home/somab/links/projects/rrg-rbeiko/somab/Proj-Thalassolituus/All_Thalassolituus_analysis/BLAST_analyses/BLASTp_NR_results"
mkdir -p "$nr_output"

for query_file in "$nr_queries"/*.faa; do
    name=$(basename "$query_file" .faa)
    blastp -query "$query_file" -db "$nr_db" -out "$nr_output/${name}_nr_blast.txt" -evalue 1e-10 -max_target_seqs 300 -num_threads "$SLURM_CPUS_PER_TASK" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs staxids sscinames sskingdoms stitle"
done

#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from pathlib import Path
from statistics import mean, median
import time

# Downloading RefSeq NifH proteins for the expanded NifH phylogeny in the manuscript, this search produced 4550 sequences at the time the analysis was performed

Entrez.email = "sm871954@dal.ca"
Entrez.tool = "Thalassolituus_nifH_analysis"

taxa = ["Bacteria", "Archaea"]
min_length = 250
max_length = 350

output_dir = Path("refseq_nifH_proteins")
output_dir.mkdir(exist_ok=True)

def fetch_proteins(taxon):
    query = f"nifH[Gene Name] AND {taxon}[Organism] AND srcdb_refseq[PROP] NOT partial[Title]"
    output_file = output_dir / f"nifH_{taxon.lower()}_refseq.faa"

    handle = Entrez.esearch(db="protein", term=query, usehistory="y", retmax=0)
    result = Entrez.read(handle)
    handle.close()

    total = int(result["Count"])
    webenv = result["WebEnv"]
    query_key = result["QueryKey"]

    print(f"Searching for NifH in {taxon}: {total} records found")

    written = 0
    batch_size = 200

    with open(output_file, "w") as output:
        for start in range(0, total, batch_size):
            handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)

            for record in SeqIO.parse(handle, "fasta"):
                if min_length <= len(record.seq) <= max_length:
                    SeqIO.write(record, output, "fasta")
                    written += 1

            handle.close()
            time.sleep(0.4)

    print(f"Kept {written} sequences in {output_file}")


def combine_files():
    output_file = output_dir / "nifH_all_refseq.faa"
    input_files = [
        output_dir / "nifH_bacteria_refseq.faa",
        output_dir / "nifH_archaea_refseq.faa"
    ]

    seen = set()
    records = []

    for input_file in input_files:
        for record in SeqIO.parse(input_file, "fasta"):
            if record.id not in seen:
                records.append(record)
                seen.add(record.id)

    SeqIO.write(records, output_file, "fasta")

    lengths = [len(record.seq) for record in records]

    print(f"{output_file}: {len(records)} sequences")
    print(f"Length: {min(lengths)}-{max(lengths)} aa")
    print(f"Mean: {mean(lengths):.1f} aa; median: {median(lengths):.1f} aa")


for taxon in taxa:
    fetch_proteins(taxon)

combine_files()

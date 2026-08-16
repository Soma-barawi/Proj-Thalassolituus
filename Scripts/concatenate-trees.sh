#!/bin/bash

mkdir -p concatenated-alignments
cd nif-alignments

# Concatenate the Nif protein alignments
AMAS.py concat -i nifH_aligned.faa nifD_aligned.faa nifK_aligned.faa -f fasta -d aa --part-format raxml \ 
  -p ../concatenated-alignments/nifHDK_partitions.txt -t ../concatenated-alignments/nifHDK_concat.faa
AMAS.py concat -i nifH_aligned.faa nifD_aligned.faa nifK_aligned.faa nifE_aligned.faa nifN_aligned.faa nifB_aligned.faa -f fasta -d aa --part-format raxml \ 
  -p ../concatenated-alignments/nifHDKENB_partitions.txt -t ../concatenated-alignments/nifHDKENB_concat.faa
cd ..

# Run IQ-TREE
iqtree3 -s concatenated-alignments/nifHDK_concat.faa -p concatenated-alignments/nifHDK_partitions.txt -m MFP+MERGE -b 1000 --alrt 1000 -T AUTO --prefix nifHDK_concat_AA
iqtree3 -s concatenated-alignments/nifHDKENB_concat.faa -p concatenated-alignments/nifHDKENB_partitions.txt -m MFP+MERGE -b 1000 --alrt 1000 -T AUTO --prefix nifHDKENB_concat_AA

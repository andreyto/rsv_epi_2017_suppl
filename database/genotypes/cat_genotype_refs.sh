#!/bin/bash
set -ex
rm -f ref.tmp.fasta
for f in orig/*.fasta; do seqret $f tmp.fasta && cat tmp.fasta >> ref.tmp.fasta; done
seqmagick convert --transcribe rna2dna ref.tmp.fasta ref.fasta
#fasta_formatter -w 0 -i ref.tmp.fasta | fasta_nucleotide_changer -d | fasta_formatter -w 80 -o ref.fasta
rm ref.tmp.fasta tmp.fasta


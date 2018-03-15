#!/bin/bash
set -ex
for subtype in A B; do
    python -m MICGENT.finders_cge annotate-cge \
        --db-cov-min 0.6 \
        --pct-ident-min 80 \
        --blast-app blastn \
        --seq ../Osmt2016_RSV_${subtype}G.fasta \
        --seq-db ../../database/genotypes/${subtype}/ref.trimmed.fasta \
        --meta-db-re '^(?P<locus>[^_]+)_(?P<cluster>[^_]+)' \
        --out-annot-csv annot_${subtype}.csv \
        --out-blast-csv blast_${subtype}.csv
done


#!/bin/bash
##Extract G gene ORF - as much as available plus some flanking sequence -
##from the OUTSMART 2017 season amplicon. Amplicon contains C-term of G,
##followed by F. Amplicons also contain the following artefacts that need to be
##removed: - 5' degenerate PCR primer (post-assembly); - flanks of XXX at both ends (some sequences).
## The XXX flanks are especially troublesome because when the collection of
## amplicons is used as BLAST database, BLAST reports match coordinates starting from the
## first non-X character rather than from the first X character (undocumented behaviour
## of BLAST). Using such shifted match coordinates will lead to incorrect extraction
## of the matching region, such as non-primer region, for example.
set -ex
#--seq-db ../../../database/reference/*_${subtype}G.fasta \
cp ../../../database/reference/KX858754_AG_forward_primer_trimmed.fasta A_trim_forward_ref.fasta
cp ../../../database/reference/KX858755_BG_forward_primer_trimmed.fasta B_trim_forward_ref.fasta
cp ../../../database/adapters/pcr_reverse_primer_rc.fasta pcr_rev_primer.fasta

for subtype in A B; do
    cp ../../../database/genotypes/${subtype}/ref.trimmed.fasta ${subtype}_genotype_ref.fasta
    ## remove flanking XXX
    python -m MICGENT.converters cat-lines --files "../${subtype}/*.fasta" \
        | seqkit replace -p '^X+' -r '' -s \
        | seqkit replace -p 'X+$' -r '' -s \
        > input_${subtype}.fasta
    ## extract part downstream of forward primer
    python -m MICGENT.finders_cge annotate-cge \
        --db-cov-min 0.6 \
        --pct-ident-min 80 \
        --blast-app blastn \
        --seq input_${subtype}.fasta \
        --seq-db ${subtype}_trim_forward_ref.fasta \
        --out-annot-csv annot_trim_forward_${subtype}.txt \
        --out-blast-csv blast_trim_forward_${subtype}.txt
    python -m MICGENT.converters extract-sequence-for-annotations  \
        --fasta-files input_${subtype}.fasta \
        --stats-file annot_trim_forward_${subtype}.json \
        --pad-left 0 \
        --pad-right 100000 \
        --omit-cluster-id \
        --keep-seq-non-matching \
        annot_trim_forward_${subtype}.txt \
        trim_forward_${subtype}.fasta
    python -m MICGENT.finders_cge annotate-cge \
        --db-cov-min 0.6 \
        --pct-ident-min 90 \
        --blast-app blastn \
        --seq trim_forward_${subtype}.fasta \
        --seq-db pcr_rev_primer.fasta \
        --db-strand plus \
        --short-reference \
        --out-annot-csv annot_pcr_rev_primer_${subtype}.txt \
        --out-blast-csv blast_pcr_rev_primer_${subtype}.txt

    #python -m pdb ~/work/micgent/python/lib/MICGENT/converters.py  extract-sequence-for-annotations \
    python -m MICGENT.converters extract-sequence-for-annotations  \
        --fasta-files trim_forward_${subtype}.fasta \
        --stats-file annot_trim_forward_${subtype}.json \
        --pad-left 100000 \
        --pad-right 0 \
        --pad-right-origin left \
        --omit-cluster-id \
        --keep-seq-non-matching \
        annot_pcr_rev_primer_${subtype}.txt \
        trim_final_${subtype}.fasta
    ## IDs must match one-to-one with the input
    seqkit seq -i -n input_${subtype}.fasta > input_${subtype}.id
    seqkit seq -i -n trim_final_${subtype}.fasta > trim_final_${subtype}.id
    cmp input_${subtype}.id trim_final_${subtype}.id
    ## *nix comm utility can be used to get the diffs after sorting IDs
    ## extract G as best match to genotyping DB C-terms plus everything
    ## upstream and a bit downstream to capture the stop codon of G
    python -m MICGENT.finders_cge annotate-cge \
        --db-cov-min 0.6 \
        --pct-ident-min 80 \
        --blast-app blastn \
        --seq trim_final_${subtype}.fasta \
        --seq-db ${subtype}_genotype_ref.fasta \
        --meta-db-re '^(?P<locus>[^_]+)_(?P<cluster>[^_]+)' \
        --out-annot-csv annot_${subtype}G.txt \
        --out-blast-csv blast_${subtype}G.txt
    python -m MICGENT.converters extract-sequence-for-annotations  \
        --fasta-files trim_final_${subtype}.fasta \
        --stats-file annot_${subtype}G.json \
        --pad-left 4000 \
        --pad-right 90 \
        annot_${subtype}G.txt \
        annot_${subtype}G.fasta
    cat annot_${subtype}G.fasta | py -x 're.sub(r"^>\S+\s+(\S+-[AB])F-",r">\1G-",x)' > annot_radar_${subtype}G.fasta
    ## IDs must match one-to-one with the input (after converting back to -F suffix)
    seqkit seq -i -n annot_radar_${subtype}G.fasta \
    | py -x 're.sub(r"^(\S+-[AB])G-(.)$",r"\1F-\2",x)' > annot_radar_${subtype}G.id
    cmp input_${subtype}.id annot_radar_${subtype}G.id
        python -m MICGENT.converters cat-lines --files "../../../database/reference/*_${subtype}G.fasta,${subtype}_genotype_ref.fasta,input_${subtype}.fasta" \
        > input_ref_${subtype}G.fasta
        mafft --maxiterate 1000 --localpair --thread -1 input_ref_${subtype}G.fasta \
        > input_ref_${subtype}G.ali.fasta
        python -m MICGENT.converters cat-lines --files "../../../database/reference/*_${subtype}G.fasta,${subtype}_genotype_ref.fasta,annot_radar_${subtype}G.fasta" \
        > annot_ref_${subtype}G.fasta
        mafft --maxiterate 1000 --localpair --thread -1 annot_ref_${subtype}G.fasta \
        > annot_ref_${subtype}G.ali.fasta
done


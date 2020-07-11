#!/usr/bin/env bash
# ==============================================================
# Tomas Bruna
# Copyright 2020, Georgia Institute of Technology, USA
#
# Run Spaln on a list of genome-protein pairs
# ==============================================================

LONG_GENE=30000
LONG_PROTEIN=15000

if [ ! "$#" -eq 3 ]; then
  echo "Usage: $0 input_batch output min_exon_score"
  exit
fi

batchFile=$1
output=$2
min_exon_score=$3

binDir="$(readlink -e $(dirname "$0"))"

export ALN_TAB="$binDir/../dependencies/spaln_table"
# Reset output
echo -n "" > $output

IFS=$'\t'

while read -r -a pair; do

  nuc=${pair[0]}
  prot=${pair[1]}

  geneLength=$(stat --printf "%s" "$nuc")
  proteinLength=$(stat --printf "%s" "$prot")
  # Estimate the maximum possible possible length of the alignment, including gaps.
  alignmentLength="$(($geneLength*2))"

  # -Q3    Algorithm runs in the fast heuristic mod
  # -pw    Report result even if alignment score is below threshold prot_id
  # -S1    Dna is in the forward orientation
  # -LS    Smith-Waterman-type local alignment. This option may prune out weakly matched terminal regions.
  # -O1    Output alignment
  # -l     Number of characters per line in alignment

  mode="-Q3"

  # Mapping mode usually consumes less memory, use it for long alignments.
  if [ $geneLength -gt $LONG_GENE ] || [ $proteinLength -gt $LONG_PROTEIN ]; then
    mode="-Q7"
  fi

  # Align and directly parse the output
  "$binDir/../dependencies/spaln" $mode -LS -pw -S1 -O1 -l $alignmentLength "$nuc" "$prot" \
    2> /dev/null | "$binDir/../dependencies/spaln_boundary_scorer" -o "${nuc}_${prot}" -w 10 \
    -s "$binDir/../dependencies/blosum62.csv" -e $min_exon_score

  cat "${nuc}_${prot}" >> $output
  rm "${nuc}_${prot}" "$nuc" "$prot"

done < "$batchFile"

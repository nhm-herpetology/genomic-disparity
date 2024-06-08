#!/bin/bash

awk '{print "mkdir "$1}' chromosome_download.config | sh

awk '!/^$/{print  > "Taxon_"NR".txt" }' chromosome_download.config

for f in *.txt; do
  [[ -f "$f" ]] || continue
  dir="${f%.*}"
  mv "$f" "$dir"
done

# 4. Execute download for each taxon, renaming the chromosomes with species names

cd Taxon_1

awk -v n=3 '{ for (i=n; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' Taxon_*.txt > search.txt

declare -a chromosomes="$(head -n 1 search.txt)" 

for i in ${chromosomes[@]}
do
   esearch -db nucleotide -query $i | efetch -format fasta > "$i".fasta
done

# 5. Rename chromosomes with species name and remove quotes

VAR1="$(awk '{print $2}' Taxon_*.txt)"

for f in *.fasta ; do mv -- "$f" "${VAR1}_$f" ; done

find . -name '*\"*' | while read f; do mv "$f" "${f//\"/}"; done

#!/bin/bash

echo "Welcome to Genomic Disparity Analysis... please ensure that your config file is formatted correctly"

echo "Extracting information from config file..."

awk '{print "mkdir "$1}' chromosome_download.config | sh

awk '!/^$/{print  > "Taxon_"NR".txt" }' chromosome_download.config

for f in *.txt; do
  [[ -f "$f" ]] || continue
  dir="${f%.*}"
  mv "$f" "$dir"
done

echo "Starting chromosome downloads..."

for dir in */; do
	(
	echo "Fetching $dir chromosomes"
	cd $dir

awk -v n=3 '{ for (i=n; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' Taxon_*.txt > search.txt

declare -a chromosomes="$(head -n 1 search.txt)" 

for i in ${chromosomes[@]}
do
   esearch -db nucleotide -query $i | efetch -format fasta > "$i".fasta
done

echo "Finished chromosome downloads... renaming files"

VAR1="$(awk '{print $2}' Taxon_*.txt)"

for f in *.fasta ; do mv -- "$f" "${VAR1}_$f" ; done

find . -name '*\"*' | while read f; do mv "$f" "${f//\"/}"; done

)done

echo "Chromosomes prepped for analysis... script finished."

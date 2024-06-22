#!/bin/bash

echo "Identifying chromosomes for extraction..."

sed 's/,/ /g' present_landmarks.csv > resultfile

head -n1 -q resultfile > new-file

awk -v n=3 '{ for (i=n; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' new-file > almost-there.txt

sed 's/"//g' almost-there.txt > chromosomes-to-find.txt

rm almost-there.txt
rm new-file
rm resultfile

mkdir chromosome_set

declare -a chromosomes="$(head -n 1 chromosomes-to-find.txt)" 

for i in ${chromosomes[@]}
do
  
find ./Taxon_** -name "mapped_$i.fasta.tsv" -exec cp "{}" ./chromosome_set \;

done

echo "Chromosome set extracted and placed in directory chromosome_set..."

echo "Renaming final chromosome files..."

cd chromosome_set/

for filename in *.fasta.tsv; do 
    [ -f "$filename" ] || continue
    mv "$filename" "${filename//mapped_/}"

done

cd ..

echo "Chromosome retrieval complete"

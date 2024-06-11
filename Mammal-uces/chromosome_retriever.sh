#!/bin/bash

echo "Identifying chromosomes for extraction..."

sed 's/,/ /g' present_landmarks.csv > resultfile

head -n1 -q resultfile > new-file

awk -v n=3 '{ for (i=n; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' new-file > chromosomes-to-find.txt

